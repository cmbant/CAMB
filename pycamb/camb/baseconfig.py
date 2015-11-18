import ctypes
import os.path as osp
import sys
import os
import platform

BASEDIR = osp.abspath(osp.dirname(__file__))
if platform.system() == "Windows":
    DLLNAME = 'cambdll.dll'
else:
    DLLNAME = 'camblib.so'


mock_load = True or os.environ.get('READTHEDOCS', None)

class Mock(object):
    def __getattr__(cls, name):
            return Mock()

if mock_load:
    CAMBL = osp.join(BASEDIR, DLLNAME)
    if not osp.isfile(CAMBL): sys.exit('camblib.so does not exist.\nPlease remove any old installation and install again.')
    camblib = ctypes.cdll.LoadLibrary(CAMBL)
else:
    camblib = Mock()

def set_filelocs():
    HighLExtrapTemplate = osp.join(BASEDIR, "HighLExtrapTemplate_lenspotentialCls.dat")
    if not osp.exists(HighLExtrapTemplate):
        HighLExtrapTemplate = osp.abspath(osp.join(BASEDIR, "../..", "HighLExtrapTemplate_lenspotentialCls.dat"))

    func = camblib.__handles_MOD_set_cls_template
    func.argtypes = [ctypes.c_char_p, ctypes.c_long]
    s = ctypes.create_string_buffer(HighLExtrapTemplate)
    func(s, ctypes.c_long(len(HighLExtrapTemplate)))


if not mock_load:
    set_filelocs()


class CAMBError(Exception):
    pass


class CAMB_Structure(ctypes.Structure):
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
