import ctypes
import os.path as osp
import sys
import platform

BASEDIR = osp.abspath(osp.dirname(__file__))
if platform.system() == "Windows":
    DLLNAME = 'cambdll.dll'
else:
    DLLNAME = 'camblib.so'

CAMBL = osp.join(BASEDIR, DLLNAME)
#CAMBL = r'C:\Work\Dist\git\camb\VisualStudio\CAMBdll\x64\Release\cambdll.dll'
#CAMBL = r'C:\Work\Dist\git\camb\VisualStudio\CAMBdll\Release\cambdll.dll'

if not osp.isfile(CAMBL): sys.exit('camblib.so does not exist.\nPlease remove any old installation and install again.')
camblib = ctypes.cdll.LoadLibrary(CAMBL)

def set_filelocs():
    HighLExtrapTemplate = osp.join(BASEDIR, "HighLExtrapTemplate_lenspotentialCls.dat")
    if not osp.exists(HighLExtrapTemplate):
        HighLExtrapTemplate = osp.abspath(osp.join(BASEDIR, "../..", "HighLExtrapTemplate_lenspotentialCls.dat"))

    func = camblib.__handles_MOD_set_cls_template
    func.argtypes = [ctypes.c_char_p, ctypes.c_long]
    s = ctypes.create_string_buffer(HighLExtrapTemplate)
    func(s, ctypes.c_long(len(HighLExtrapTemplate)))

set_filelocs()

class CAMBError(Exception):
    pass


class CAMB_Structure(ctypes.Structure):
    def __str__(self):
        s = ''
        for field_name, field_type in self._fields_:
            obj = getattr(self, field_name)
            if isinstance(obj, CAMB_Structure):
                s += field_name + ':\n  ' + str(obj).replace('\n', '\n  ')
            else:
                if isinstance(obj, ctypes.Array):
                    s += field_name + ' = ' + str(obj[:min(7, len(obj))]) + '\n'
                else:
                    s += field_name + ' = ' + str(obj) + '\n'
        return s
