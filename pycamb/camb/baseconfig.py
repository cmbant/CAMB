import os.path as osp
import sys
import os
import six
import platform

BASEDIR = osp.abspath(osp.dirname(__file__))
if platform.system() == "Windows":
    DLLNAME = 'cambdll.dll'

else:
    DLLNAME = 'camblib.so'
CAMBL = osp.join(BASEDIR, DLLNAME)

mock_load = os.environ.get('READTHEDOCS', None)

if not mock_load:
    import ctypes
    from ctypes import Structure, POINTER, byref


    class ifort_gfortran_loader(ctypes.CDLL):

        def __getitem__(self, name_or_ordinal):
            try:
                res = super(ifort_gfortran_loader, self).__getitem__(name_or_ordinal)
            except:
                # ifort style exports instead
                res = super(ifort_gfortran_loader, self).__getitem__(
                    name_or_ordinal.replace('_MOD_', '_mp_').replace('__', '') + '_')
            return res


    if not osp.isfile(CAMBL):
        if platform.system() == "Windows":
            # allow local git loading if not installed
            import struct

            is32Bit = struct.calcsize("P") == 4
            CAMBL = osp.join(BASEDIR, '..', 'dlls', ('cambdll_x64.dll', DLLNAME)[is32Bit])
        if not osp.isfile(CAMBL):
            sys.exit('%s does not exist.\nPlease remove any old installation and install again.' % DLLNAME)
    camblib = ctypes.LibraryLoader(ifort_gfortran_loader).LoadLibrary(CAMBL)
else:
    # This is just so readthedocs build will work without CAMB binary library
    try:
        from unittest.mock import MagicMock
    except ImportError:
        from mock import Mock as MagicMock


    class Mock(MagicMock):
        @classmethod
        def __getattr__(cls, name):
            if name == 'pi':
                return 1
            else:
                return Mock()

        def __mul__(self, other):
            return Mock()

        def __pow__(self, other):
            return 1


    MOCK_MODULES = ['numpy', 'numpy.ctypeslib', 'ctypes']
    sys.modules.update((mod_name, Mock()) for mod_name in MOCK_MODULES)
    camblib = Mock()
    Structure = object
    import ctypes


def dll_import(tp, module, func):
    try:
        # gfortran
        return tp.in_dll(camblib, "__%s_MOD_%s" % (module, func))
    except:
        # ifort
        return tp.in_dll(camblib, "%s_mp_%s_" % (module, func))


def set_filelocs():
    HighLExtrapTemplate = osp.join(BASEDIR, "HighLExtrapTemplate_lenspotentialCls.dat")
    if not osp.exists(HighLExtrapTemplate):
        HighLExtrapTemplate = osp.abspath(osp.join(BASEDIR, "../..", "HighLExtrapTemplate_lenspotentialCls.dat"))
    HighLExtrapTemplate = six.b(HighLExtrapTemplate)
    func = camblib.__handles_MOD_set_cls_template
    func.argtypes = [ctypes.c_char_p, ctypes.c_long]
    s = ctypes.create_string_buffer(HighLExtrapTemplate)
    func(s, ctypes.c_long(len(HighLExtrapTemplate)))


if not mock_load:
    set_filelocs()
    from ctypes import c_byte, c_int, c_bool

    _get_allocatable_size = camblib.__handles_MOD_getallocatablesize
    _get_allocatable_size.restype = c_int
    f_allocatable = c_byte * (_get_allocatable_size() // 8)
else:
    f_allocatable = None


class CAMBError(Exception):
    pass


class CAMBValueError(ValueError):
    pass


class CAMBUnknownArgumentError(ValueError):
    pass


class CAMBParamRangeError(CAMBError):
    pass


class CAMB_Structure(Structure):
    def __str__(self):
        s = ''
        for field_name, field_type in self._fields_:
            obj = getattr(self, field_name)
            if field_name[0:2] == '__': continue
            if field_name[0] == '_':
                field_name = field_name[1:]
                obj = getattr(self, field_name)
            if isinstance(obj, CAMB_Structure):
                s += field_name + ':\n  ' + str(obj).replace('\n', '\n  ').strip(' ')
            else:
                if isinstance(obj, ctypes.Array):
                    s += field_name + ' = ' + str(obj[:min(7, len(obj))]) + '\n'
                else:
                    s += field_name + ' = ' + str(obj) + '\n'
        return s


class F2003Class(CAMB_Structure):
    # Wraps a fortran type (potentially containing allocatable elements that can't be accessed directly by ctypes)
    # Assumes fortran subroutines (handles module in camb_python) called CLASS_new and CLASS_free exist which actually
    # allocate and deallocate instances the corresponding fortran types.

    _instance_count = {}
    _imports = {}

    def __new__(cls, *args, **kwargs):
        return cls._new_copy()

    @classmethod
    def _new_copy(cls, source=None):
        if source is None:
            _key = POINTER(cls)()
        else:
            _key = POINTER(cls)(source)
        new = cls.import_func('new', pointer=True)
        new(byref(_key))
        instance = _key.contents
        instance._key = _key
        cls._instance_count[cls] = cls._instance_count.get(cls, 0) + 1
        return instance

    def copy(self):
        """
        Make independent copy of this object.

        :return: copy of self
        """
        res = self._new_copy(source=self)
        res._init_members()
        return res

    def __init__(self, **kwargs):
        self._init_members(**kwargs)

    @classmethod
    def import_func(cls, tag, pointer=False, extra_args=[], restype=None):
        # Import fortran function called CLASSNAME_tag
        func = cls._imports.get((cls, tag), None)
        if func is None:
            func = getattr(camblib, '__handles_MOD_' + cls.__name__.lower() + '_' + tag)
            if pointer:
                func.argtypes = [POINTER(POINTER(cls))] + extra_args
            else:
                func.argtypes = [POINTER(cls)] + extra_args
            if restype is not None: func.restype = restype
            cls._imports[(cls, tag)] = func
        return func

    def call_func(self, tag, extra_args=[], args=[], restype=None):
        func = self.import_func(tag, extra_args=extra_args, restype=restype)
        return func(byref(self), *args)

    def _init_members(self, **kwargs):
        # e.g. set up allocatables or other things not directly accessible
        # called both on new instance and on making a copy
        pass

    def free_instance(self):
        cls = self.__class__
        free = cls.import_func('free', pointer=True)
        free(byref(self._key))
        cls._instance_count[cls] -= 1
        if cls._instance_count[cls] == 0:
            del cls._instance_count[cls]

    def __del__(self):
        if hasattr(self, '_key'):
            self.free_instance()

    def _void_p(self):
        return ctypes.cast(ctypes.pointer(self), ctypes.c_void_p)
