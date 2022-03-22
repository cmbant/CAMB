import os
import os.path as osp
import sys
import platform
import ctypes
from ctypes import Structure, POINTER, byref, c_int, c_double, c_bool, c_float
from numpy.ctypeslib import ndpointer
import numpy as np


def ndpointer_or_null(*args, **kwargs):
    # allows passing None to fortran optional arguments
    # from https://stackoverflow.com/a/37664693/1022775
    base = ndpointer(*args, **kwargs)

    def from_param(cls, obj):
        if obj is None:
            return obj
        return base.from_param(obj)

    return type(base.__name__, (base,), {'from_param': classmethod(from_param)})


numpy_3d = ndpointer(c_double, flags='C_CONTIGUOUS', ndim=3)
numpy_2d = ndpointer(c_double, flags='C_CONTIGUOUS', ndim=2)
numpy_1d = ndpointer(c_double, flags='C_CONTIGUOUS')
numpy_1d_or_null = ndpointer_or_null(c_double, flags='C_CONTIGUOUS')
numpy_1d_int = ndpointer(c_int, flags='C_CONTIGUOUS')

BASEDIR = osp.abspath(osp.dirname(__file__))
if platform.system() == "Windows":
    DLLNAME = 'cambdll.dll'

else:
    DLLNAME = 'camblib.so'
CAMBL = osp.join(BASEDIR, DLLNAME)

gfortran = True


class IfortGfortranLoader(ctypes.CDLL):

    def __getitem__(self, name_or_ordinal):
        if gfortran:
            res = super().__getitem__(name_or_ordinal)
        else:
            res = super().__getitem__(name_or_ordinal.replace('_MOD_', '_mp_').replace('__', '') + '_')
        return res


mock_load = os.environ.get('CAMB_MOCK_LOAD', None)
if mock_load:
    # noinspection PyCompatibility
    from unittest.mock import MagicMock

    camblib = MagicMock()
    import_property = MagicMock()
else:

    if not osp.isfile(CAMBL):
        sys.exit(
            'Library file %s does not exist.\nMake sure you have installed or built the camb package '
            '(e.g. using "python setup.py make"); or remove any old conflicting installation and install again.'
            % CAMBL)

    camblib = ctypes.LibraryLoader(IfortGfortranLoader).LoadLibrary(CAMBL)

    try:
        c_int.in_dll(camblib, "handles_mp_set_cls_template_")
        gfortran = False
    except Exception:
        pass


    class _dll_value:
        __slots__ = ['f']

        def __init__(self, f):
            self.f = f

        def __get__(self, instance, owner):
            return self.f.value

        def __set__(self, instance, value):
            self.f.value = value


    def import_property(tp, module, func):
        if gfortran:
            f = tp.in_dll(camblib, "__%s_MOD_%s" % (module.lower(), func.lower()))
        else:
            f = tp.in_dll(camblib, "%s_mp_%s_" % (module.lower(), func.lower()))

        return _dll_value(f)


def lib_import(module_name, class_name, func_name, restype=None):
    if class_name:
        class_name += '_'
    func = getattr(camblib, '__' + module_name.lower() +
                   '_MOD_' + (class_name + func_name).lower())
    if restype:
        func.restype = restype
    return func


def set_cl_template_file(cl_template_file=None):
    if cl_template_file and not osp.exists(cl_template_file):
        raise ValueError('File not found : %s' % cl_template_file)

    template = cl_template_file or osp.join(BASEDIR,
                                            "HighLExtrapTemplate_lenspotentialCls.dat")
    if not osp.exists(template):
        template = osp.abspath(
            osp.join(BASEDIR, "..", "fortran",
                     "HighLExtrapTemplate_lenspotentialCls.dat"))
    template = template.encode("latin-1")
    func = camblib.__handles_MOD_set_cls_template
    func.argtypes = [ctypes.c_char_p, ctypes.c_long]
    s = ctypes.create_string_buffer(template)
    func(s, ctypes.c_long(len(template)))


def check_fortran_version(version):
    if mock_load:
        return
    func = camblib.__camb_MOD_camb_getversion
    func.argtypes = [ctypes.c_char_p, ctypes.c_long]
    s = ctypes.create_string_buffer(33)
    func(s, ctypes.c_long(32))
    fortran_version = s.value.decode('ascii').strip()
    if fortran_version != version:
        raise CAMBFortranError(
            'Version %s of fortran library does not match python version (%s).' %
            (fortran_version,
             version) + '\nUpdate install or use "setup.py make" to rebuild library.'
            + '\n(also check camb.__file__ is actually at the path you think you are loading)')


set_cl_template_file()


def _get_fortran_sizes():
    _get_allocatable_size = camblib.__handles_MOD_getallocatablesize
    allocatable = c_int()
    allocatable_array = c_int()
    allocatable_object_array = c_int()
    _get_allocatable_size(byref(allocatable), byref(allocatable_array),
                          byref(allocatable_object_array))
    return allocatable.value, allocatable_array.value, allocatable_object_array.value


_f_allocatable_size, _f_allocatable_array_size, _f_allocatable_object_array_size = _get_fortran_sizes()
assert _f_allocatable_size % ctypes.sizeof(ctypes.c_void_p) == 0 and \
       _f_allocatable_array_size % ctypes.sizeof(ctypes.c_void_p) == 0 and \
       _f_allocatable_object_array_size % ctypes.sizeof(ctypes.c_void_p) == 0

# make dummy type of right size to hold fortran allocatable; must be ctypes pointer type to keep auto-alignment correct
f_pointer = ctypes.c_void_p * (_f_allocatable_size // ctypes.sizeof(ctypes.c_void_p))

# These are used for general types, so avoid type checking (and making problematic temporary cast objects)
_get_allocatable = lib_import('handles', 'F2003Class', 'GetAllocatable')
_set_allocatable = lib_import('handles', 'F2003Class', 'SetAllocatable')
_new_instance = lib_import('handles', 'F2003Class', 'new')
_free_instance = lib_import('handles', 'F2003Class', 'free')
_get_id = lib_import('handles', 'F2003Class', 'get_id')


class FortranAllocatable(Structure):
    pass


_reuse_pointer = ctypes.c_void_p()
_reuse_typed_id = f_pointer()


# member corresponding to class(...), allocatable :: member in fortran
class _AllocatableObject(FortranAllocatable):
    _fields_ = [("allocatable", f_pointer)]

    def get_allocatable(self):
        _get_allocatable(byref(self), byref(_reuse_typed_id), byref(_reuse_pointer))
        if _reuse_pointer:
            return ctypes.cast(_reuse_pointer, POINTER(
                F2003Class._class_pointers[tuple(_reuse_typed_id)])).contents
        else:
            return None

    def set_allocatable(self, instance, name):
        if instance and not isinstance(instance, self._baseclass):
            raise TypeError(
                '%s expects object that is an instance of %s' % (name, self._baseclass.__name__))
        _set_allocatable(byref(self), byref(instance.fortran_self) if instance else None)


_class_cache = {}


# noinspection PyPep8Naming
def AllocatableObject(cls=None):
    if cls is None:
        cls = F2003Class
    if not issubclass(cls, F2003Class):
        raise ValueError("AllocatableObject type must be descended from F2003Class")
    res = _class_cache.get(cls, None)
    if res:
        return res
    else:
        res = type("Allocatable" + cls.__name__, (_AllocatableObject,),
                   {"_baseclass": cls})
        _class_cache[cls] = res
        return res


class _AllocatableArray(FortranAllocatable):  # member corresponding to allocatable :: d(:) member in fortran
    _fields_ = [("allocatable", ctypes.c_void_p * (
            _f_allocatable_array_size // ctypes.sizeof(ctypes.c_void_p)))]

    def get_allocatable(self):
        size = self._get_allocatable_1D_array(byref(self), byref(_reuse_pointer))
        if size:
            return ctypes.cast(_reuse_pointer, POINTER(self._ctype * size)).contents
        else:
            return np.empty(0)

    def set_allocatable(self, array, name):
        self._set_allocatable_1D_array(byref(self), np.array(array, dtype=self._dtype),
                                       byref(c_int(0 if array is None else len(array))))


class _ArrayOfAllocatable(FortranAllocatable):
    def __getitem__(self, item):
        value = self.allocatables[item]
        if isinstance(value, list):
            return [x.get_allocatable() for x in value]
        else:
            return value.get_allocatable()

    def __setitem__(self, key, value):
        alloc = self.allocatables[key]
        alloc.set_allocatable(value, self.__class__.__name__)

    def __len__(self):
        return len(self.allocatables)

    def __repr__(self):
        s = ''
        for i in range(len(self.allocatables)):
            item = self[i]
            content = item._as_string() if isinstance(item, CAMB_Structure) else repr(item)
            s += ('%s: <%s>\n  ' % (i, item.__class__.__name__) + content.replace('\n', '\n  ')).strip(' ')
        return s


def _make_array_class(baseclass, size):
    res = _class_cache.get((baseclass, size), None)
    if res:
        return res

    class Temp(_ArrayOfAllocatable):
        _fields_ = [("allocatables", AllocatableObject(baseclass) * size)]

    Temp.__name__ = "%sArray_%s" % (baseclass.__name__, size)
    _class_cache[(baseclass, size)] = Temp
    return Temp


class _AllocatableObjectArray(FortranAllocatable):
    # member corresponding to allocatable :: d(:) array of allocatable classes
    _fields_ = [("allocatable", ctypes.c_void_p * (
            _f_allocatable_object_array_size // ctypes.sizeof(ctypes.c_void_p)))]

    def get_allocatable(self):
        size = self._get_allocatable_object_1D_array(byref(self), byref(_reuse_pointer))
        if size:
            return ctypes.cast(_reuse_pointer,
                               POINTER(_make_array_class(self._baseclass, size))).contents
        else:
            return []

    def set_allocatable(self, array, name):
        if array is None:
            array = []
        pointers = (f_pointer * len(array))()
        for i, instance in enumerate(array):
            if not isinstance(instance, self._baseclass):
                raise TypeError(
                    '%s expects object that is an instance of %s' % (name, self._baseclass.__name__))
            pointers[i] = instance.fortran_self

        self._set_allocatable_object_1D_array(byref(self), byref(pointers), byref(c_int(len(array))))


_AllocatableObjectArray._get_allocatable_object_1D_array = camblib.__handles_MOD_get_allocatable_object_1d_array
_AllocatableObjectArray._get_allocatable_object_1D_array.restype = c_int
_AllocatableObjectArray._set_allocatable_object_1D_array = camblib.__handles_MOD_set_allocatable_object_1d_array


# noinspection PyPep8Naming
def AllocatableObjectArray(cls=None):
    if cls is None:
        cls = F2003Class
    if not issubclass(cls, F2003Class):
        raise ValueError("AllocatableObject type must be descended from F2003Class")
    return type("AllocatableArray" + cls.__name__, (_AllocatableObjectArray,),
                {"_baseclass": cls})


class AllocatableArrayInt(_AllocatableArray):
    _dtype = int
    _ctype = c_int


AllocatableArrayInt._get_allocatable_1D_array = camblib.__handles_MOD_get_allocatable_1d_array_int
AllocatableArrayInt._get_allocatable_1D_array.restype = c_int
AllocatableArrayInt._set_allocatable_1D_array = camblib.__handles_MOD_set_allocatable_1d_array_int
AllocatableArrayInt._set_allocatable_1D_array.argtypes = [POINTER(AllocatableArrayInt), numpy_1d_int, POINTER(c_int)]


class AllocatableArrayDouble(_AllocatableArray):
    _dtype = np.float64
    _ctype = c_double


AllocatableArrayDouble._get_allocatable_1D_array = camblib.__handles_MOD_get_allocatable_1d_array
AllocatableArrayDouble._get_allocatable_1D_array.restype = c_int
AllocatableArrayDouble._set_allocatable_1D_array = camblib.__handles_MOD_set_allocatable_1d_array
AllocatableArrayDouble._set_allocatable_1D_array.argtypes = [POINTER(AllocatableArrayDouble), numpy_1d, POINTER(c_int)]


def fortran_array(c_pointer, shape, dtype=np.float64, order='F', own_data=True):
    if not hasattr(shape, '__len__'):
        shape = np.atleast_1d(shape)
    arr_size = np.prod(shape[:]) * np.dtype(dtype).itemsize
    buf_from_mem = ctypes.pythonapi.PyMemoryView_FromMemory
    buf_from_mem.restype = ctypes.py_object
    buf_from_mem.argtypes = (ctypes.c_void_p, ctypes.c_int, ctypes.c_int)
    buffer = buf_from_mem(c_pointer, arr_size, 0x100)
    arr = np.ndarray(tuple(shape[:]), dtype, buffer, order=order)
    if own_data and not arr.flags.owndata:
        return arr.copy()
    else:
        return arr


class CAMBError(Exception):
    pass


class CAMBValueError(ValueError):
    pass


class CAMBUnknownArgumentError(ValueError):
    pass


class CAMBParamRangeError(CAMBError):
    pass


class CAMBFortranError(Exception):
    pass


def method_import(module_name, class_name, func_name, restype=None, extra_args=(),
                  nopass=False):
    func = lib_import(module_name, class_name, func_name, restype)
    if extra_args is not None and len(extra_args):
        func.argtypes = ([] if nopass else [POINTER(f_pointer)]) + list(extra_args)
    return func


# Handle custom field types inspired by:
# https://stackoverflow.com/questions/45527945/extend-ctypes-to-specify-field-overloading
class FortranManagedField:
    __slots__ = ['name', 'real_name', 'type_']

    def __init__(self, name, type_):
        self.name = name
        self.type_ = type_
        self.real_name = "_" + name

    def __get__(self, instance, owner):
        value = getattr(instance, self.real_name)
        if issubclass(self.type_, FortranAllocatable):
            return value.get_allocatable()
        return value

    def __set__(self, instance, value):
        field_value = getattr(instance, self.real_name)
        if issubclass(self.type_, FortranAllocatable):
            field_value.set_allocatable(value, self.name)
            return
        if issubclass(self.type_, F2003Class):
            field_value.replace(value)
            return
        setattr(instance, self.real_name, value)


class NamedIntField:
    __slots__ = ['real_name', 'values', 'name_values']

    def __init__(self, name, **kwargs):
        self.real_name = "_" + name
        names = kwargs["names"]
        self.values = {}
        if isinstance(names, (list, tuple)):
            self.name_values = {}
            start = kwargs.get("start", 0)
            for i, name in enumerate(names):
                self.name_values[name] = i + start
                self.values[i + start] = name
        else:
            assert isinstance(names, dict)
            self.name_values = names
            for name in names:
                self.values[names[name]] = name

    def __get__(self, instance, owner):
        value = getattr(instance, self.real_name)
        return self.values[value]

    def __set__(self, instance, value):
        if isinstance(value, str):
            value = self.name_values[value]
        elif value not in self.values:
            raise ValueError("Value %s not in allowed: %s" % (value, self.name_values))
        setattr(instance, self.real_name, value)


class BoolField:  # fortran-compatible boolean (actually c_int internally)
    __slots__ = ['real_name']

    def __init__(self, name):
        self.real_name = "_" + name

    def __get__(self, instance, owner):
        return getattr(instance, self.real_name) != 0

    def __set__(self, instance, value):
        setattr(instance, self.real_name, (0, 1)[value])


class SizedArrayField:  # statically sized array with another field determining size
    __slots__ = ['real_name', 'size_name']

    def __init__(self, name, size_name):
        self.real_name = "_" + name
        self.size_name = size_name

    def __get__(self, instance, owner):
        size = getattr(instance, self.size_name)
        value = getattr(instance, self.real_name)
        if size == len(value):
            return value
        return POINTER(value._type_ * size)(value).contents

    def __set__(self, instance, value):
        field = getattr(instance, self.real_name)
        if len(value) > len(field):
            raise CAMBParamRangeError(
                "%s can be of max length %s" % (self.real_name[1:], len(field)))
        field[:len(value)] = value
        setattr(instance, self.size_name, len(value))


class CAMBStructureMeta(type(Structure)):
    # noinspection PyMethodParameters
    def __new__(metacls, name, bases, namespace):
        _fields = namespace.get("_fields_", "")
        ctypes_fields = []
        try:
            F2003 = F2003Class
        except NameError:
            class F2003:
                pass

        tps = {c_bool: "boolean", c_double: "float64", c_int: "integer",
               c_float: "float32",
               AllocatableArrayDouble: "float64 array",
               AllocatableArrayInt: "integer array",
               ctypes.c_void_p: "pointer"}
        field_doc = ''
        for field in _fields:
            field_name = field[0]
            field_type = field[1]
            if field_type == c_bool:
                new_field = BoolField(field_name)
                ctypes_fields.append(("_" + field_name, c_int))
            elif issubclass(field_type, FortranAllocatable) or issubclass(field_type, F2003):
                new_field = FortranManagedField(field_name, field_type)
                ctypes_fields.append(("_" + field_name, field_type))
            elif len(field) > 2 and isinstance(field[2], dict):
                dic = field[2]
                if "names" in dic:
                    if field_type != c_int:
                        raise CAMBFortranError("Named fields only allowed for c_int")
                    new_field = NamedIntField(field_name, **dic)
                    ctypes_fields.append(("_" + field_name, c_int))
                elif "size" in dic:
                    if not issubclass(field_type, ctypes.Array):
                        raise CAMBFortranError(
                            "sized fields only allowed for ctypes Arrays")
                    if dic["size"] not in [x[0] for x in _fields]:
                        raise CAMBFortranError(
                            "size must be name of field in same structure (%s for %s)" % (
                                dic["size"], field_name))
                    new_field = SizedArrayField(field_name, dic["size"])
                    ctypes_fields.append(("_" + field_name, field_type))
                else:
                    raise CAMBFortranError(
                        "Unknown dictionary content for %s, %s" % (field_name, dic))
            else:
                new_field = None
            if new_field:
                namespace[field_name] = new_field
            else:
                ctypes_fields.append((field_name, field_type))
            if field[0][0] != '_':  # add :ivar: documentation for each field
                field_doc += "\n    :ivar %s:" % field[0]
                if isinstance(field[-1], dict) and field[-1].get('names', None):
                    field_doc += " (integer/string, one of: %s) " % (", ".join(field[-1]["names"]))
                else:
                    tp = tps.get(field[1], None)
                    if tp:
                        field_doc += " (*%s*)" % tp
                    elif issubclass(field[1], ctypes.Array):
                        field_doc += " (*%s array*)" % tps[field[1]._type_]
                    elif issubclass(field[1], CAMB_Structure):
                        field_doc += " :class:`%s.%s`" % (field[1].__module__, field[1].__name__)
                    elif issubclass(field[1], _AllocatableObject):
                        field_doc += " :class:`%s.%s`" % (
                            field[1]._baseclass.__module__, field[1]._baseclass.__name__)
                    elif issubclass(field[1], _AllocatableObjectArray):
                        field_doc += " array of :class:`%s.%s`" % (
                            field[1]._baseclass.__module__, field[1]._baseclass.__name__)

                if len(field) > 2 and not isinstance(field[-1], dict):
                    field_doc += " " + field[-1]

        namespace["_fields_"] = ctypes_fields

        if field_doc:
            namespace['__doc__'] = namespace.get('__doc__', "") + "\n" + field_doc

        # noinspection PyTypeChecker
        cls: CAMB_Structure = super().__new__(metacls, name, bases, namespace)

        if name == "F2003Class" or issubclass(bases[0], F2003):
            cls._class_imports = {}
            if "_fortran_class_name_" not in cls.__dict__:
                cls._fortran_class_name_ = name

        prefix = getattr(cls, '_method_prefix_', "f_")
        methods = cls.__dict__.get("_methods_", "")

        def make_method(_func, _name, _nopass, doc):
            if _nopass:
                def method_func(self, *args):
                    return _func(*args)
            else:
                def method_func(self, *args):
                    return _func(self.fortran_self, *args)
            if doc:
                method_func.__doc__ = doc
            method_func.__name__ = _name
            return method_func

        for method in methods:
            method_name = method[0]
            extra_args = method[1]
            restype = method[2] if len(method) > 2 else None
            opts = method[3] if len(method) > 3 else {}
            nopass = opts.get("nopass", False)
            try:
                func = method_import(cls._fortran_class_module_, cls._fortran_class_name_,
                                     method_name, extra_args=extra_args, nopass=nopass, restype=restype)
            except AttributeError:
                raise AttributeError('No function %s_%s found in module %s' % (
                    cls._fortran_class_name_, method_name, cls._fortran_class_module_))
            new_method = make_method(func, prefix + method_name, nopass, opts.get("doc", ""))
            setattr(cls, prefix + method_name, new_method)

        return cls


# noinspection PyPep8Naming
class CAMB_Structure(Structure, metaclass=CAMBStructureMeta):

    # noinspection PyUnresolvedReferences
    @classmethod
    def get_all_fields(cls):
        if cls != CAMB_Structure:
            fields = cls.__bases__[0].get_all_fields()
        else:
            fields = []
        fields += [(name[1:], value) if name.startswith('_') else (name, value) for name, value in
                   cls.__dict__.get('_fields_', []) if
                   not name.startswith('__')]
        return fields

    @classmethod
    def get_valid_field_names(cls):
        return set(field[0] for field in cls.get_all_fields())

    def _as_string(self):
        s = ''
        for field_name, field_type in self.get_all_fields():
            obj = getattr(self, field_name)
            if isinstance(obj, (CAMB_Structure, FortranAllocatable)):
                content = obj._as_string() if isinstance(obj, CAMB_Structure) else str(obj)
                s += (field_name + ': <%s>\n  ' % obj.__class__.__name__ + content.replace('\n', '\n  ')).strip(' ')
            else:
                if isinstance(obj, ctypes.Array):
                    if len(obj) > 20:
                        s += field_name + ' = ' + str(obj[:7])[:-1] + ', ...]\n'
                    else:
                        s += field_name + ' = ' + str(obj[:len(obj)]) + '\n'
                else:
                    s += field_name + ' = ' + str(obj) + '\n'
        return s

    def __repr__(self):
        return 'class: <%s>\n ' % self.__class__.__name__ + self._as_string().replace('\n', '\n ')


class _FortranSelf:
    def __get__(self, instance, owner):
        if not instance and owner.__class__ is CAMBStructureMeta:
            # prevent error during introspection of classes
            return None
        pointer = f_pointer()
        owner._fortran_selfpointer_function(byref(ctypes.pointer(instance)), byref(pointer))
        instance.fortran_self = pointer
        return pointer


class F2003Class(CAMB_Structure):
    # Wraps a fortran type/class that is allocated in fortran, potentially containing allocatable  _fields_
    # elements that are instances of classes, allocatable arrays that are wrapped in python, and list
    # of class _methods_ from fortran.
    #
    # Note that assigning to allocatable fields makes a deep copy of the object so the object always owns all memory
    # belonging to its fields. Accessing an allocatable field makes a new class pointer object on the fly. It can
    # become undefined if the allocatable field is reassigned.

    # classes are referenced by their fortran null pointer object. _class_pointers is a dictionary relating these
    # f_pointer to python classes. Elements are added each class by the @fortran_class decorator.
    _class_pointers = {}

    # dictionary mapping class names to classes
    _class_names = {}

    __slots__ = ()

    # pointer to fortran class; generated once per instance using _fortran_selfpointer_function then replaced with
    # the actual value
    fortran_self = _FortranSelf()

    def __new__(cls, *args, **kwargs):
        return cls._new_copy()

    def __init__(self, **kwargs):
        unknowns = set(kwargs) - self.get_valid_field_names()
        if unknowns:
            raise ValueError('Unknown argument(s): %s' % unknowns)
        super().__init__(**kwargs)

    @classmethod
    def _new_copy(cls, source=None):
        if source is None:
            _key = POINTER(cls)()
        else:
            _key = POINTER(cls)(source)
        pointer_func = getattr(cls, '_fortran_selfpointer_function', None)
        if pointer_func is None:
            if getattr(cls, '_optional_compile', False):
                raise CAMBFortranError(
                    'Class %s has not been built into the Fortran binary,'
                    ' edit Makefile_main and rebuild to use.' % cls.__name__)
            raise CAMBFortranError(
                'Cannot instantiate %s, is base class or needs @fortran_class decorator' % cls.__name__)
        _new_instance(byref(_key), byref(pointer_func))
        instance = _key.contents
        instance._key = _key
        return instance

    def copy(self):
        """
        Make independent copy of this object.

        :return: deep copy of self
        """
        return self._new_copy(source=self)

    __copy__ = copy

    @classmethod
    def import_method(cls, tag, extra_args=(), restype=None, nopass=False,
                      allow_inherit=True):
        func = cls._class_imports.get(tag, None)
        if func is None:
            try:
                func = method_import(cls._fortran_class_module_, cls._fortran_class_name_,
                                     tag, extra_args=extra_args, nopass=nopass, restype=restype)
            except AttributeError:
                try:
                    if not allow_inherit or cls.__bases__[0] == F2003Class:
                        raise
                    # noinspection PyUnresolvedReferences
                    func = cls.__bases__[0].import_method(tag, extra_args, restype, nopass=nopass)
                except AttributeError:
                    raise AttributeError(
                        'No function %s_%s found ' % (cls._fortran_class_name_, tag))
            cls._class_imports[tag] = func
        return func

    def call_method(self, tag, extra_args=(), args=(), restype=None, nopass=False, allow_inherit=True):
        func = self.import_method(tag, extra_args=extra_args, restype=restype,
                                  nopass=nopass, allow_inherit=allow_inherit)
        if nopass:
            return func(*args)
        else:
            return func(byref(self.fortran_self), *args)

    def __del__(self):
        key = getattr(self, '_key', None)
        if key:
            _free_instance(byref(key), byref(self.__class__._fortran_selfpointer_function))

    def replace(self, instance):
        """
        Replace the content of this class with another instance, doing a deep copy (in Fortran)

        :param instance: instance of the same class to replace this instance with
        """
        if type(instance) != type(self):
            raise TypeError(
                'Cannot assign non-identical types (%s to %s, non-allocatable)' % (type(instance), type(self)))
        self.call_method('Replace', extra_args=[POINTER(f_pointer)], allow_inherit=False,
                         args=[byref(instance.fortran_self)])

    @staticmethod
    def make_class_named(name, base_class=None):
        if not isinstance(name, type):
            cls = F2003Class._class_names.get(name, None)
            if not cls:
                raise CAMBValueError("Class not found: %s" % name)
        else:
            cls = name
        if base_class is None or issubclass(cls, base_class):
            return cls()
        else:
            raise CAMBValueError(
                "class %s is not a type of %s" % (cls.__name__, base_class.__name__))


# Decorator to get function to get class pointers to each class type, and build index of classes
# that allocatables could have
def fortran_class(cls, optional=False):
    if mock_load:
        return cls
    class_module = getattr(cls, "_fortran_class_module_", None)
    if not class_module:
        msg = "F2003Class %s must define _fortran_class_module_" % cls.__name__
        print(msg)
        raise CAMBFortranError(msg)
    try:
        cls._fortran_selfpointer_function = lib_import(class_module, cls._fortran_class_name_,
                                                       'selfpointer')
    except AttributeError as e:
        if optional:
            cls._optional_compile = True
            return cls
        else:
            print(e)
            raise CAMBFortranError("Class %s cannot find fortran %s_SelfPointer method in module %s." %
                                   (cls.__name__, cls._fortran_class_name_, class_module))

    _get_id(byref(cls._fortran_selfpointer_function), byref(_reuse_typed_id))
    F2003Class._class_pointers[tuple(_reuse_typed_id)] = cls
    F2003Class._class_names[cls.__name__] = cls
    return cls


def optional_fortran_class(cls):
    return fortran_class(cls, optional=True)
