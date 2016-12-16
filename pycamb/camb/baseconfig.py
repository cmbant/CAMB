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
    from ctypes import Structure


    class ifort_gfortran_loader(ctypes.CDLL):

        def __getitem__(self, name_or_ordinal):
            try:
                res = super(ifort_gfortran_loader, self).__getitem__(name_or_ordinal)
            except:
                # ifort style exports instead
                res = super(ifort_gfortran_loader, self).__getitem__(
                    name_or_ordinal.replace('_MOD_', '_mp_').replace('__', '') + '_')
            return res


    if not osp.isfile(CAMBL): sys.exit(
        '%s does not exist.\nPlease remove any old installation and install again.' % DLLNAME)
    camblib = ctypes.LibraryLoader(ifort_gfortran_loader).LoadLibrary(CAMBL)
# camblib = ctypes.cdll.LoadLibrary(CAMBL)
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


class CAMBError(Exception):
    pass


class CAMB_Structure(Structure):
    def __str__(self):
        s = ''
        for field_name, field_type in self._fields_:
            obj = getattr(self, field_name)
            if isinstance(obj, CAMB_Structure):
                s += field_name + ':\n  ' + str(obj).replace('\n', '\n  ').strip(' ')
            else:
                if isinstance(obj, ctypes.Array):
                    s += field_name + ' = ' + str(obj[:min(7, len(obj))]) + '\n'
                else:
                    s += field_name + ' = ' + str(obj) + '\n'
        return s
