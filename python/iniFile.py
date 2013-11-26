# AL Apr 11,Jun13
import os

class iniFile:


    def __init__(self, filename=None, keep_includes=False):

        self.params = dict()
        self.readOrder = []
        self.defaults = []
        self.includes = []
        if filename is not None and filename!='': self.readFile(filename, keep_includes)


    def readFile(self, filename, keep_includes=False, if_not_defined=False):
        fileincludes = []
        filedefaults = []
        textFileHandle = open(filename)
        # Remove blanck lines and comment lines from the python list of lists.
        for line in textFileHandle:
            s = line.strip()
            if s == 'END':break
            if s.startswith('#'):continue
            elif s.startswith('INCLUDE('):
                fileincludes.append(s[s.find('(') + 1:s.rfind(')')])
            elif s.startswith('DEFAULT('):
                filedefaults.append(s[s.find('(') + 1:s.rfind(')')])
            elif s != '':
                eq = s.find('=')
                if eq >= 0:
                    key = s[0:eq].strip()
                    if key in self.params:
                        if if_not_defined: continue
                        raise Exception('Error: duplicate key: ' + key)
                    value = s[eq + 1:].strip()
                    self.params[key] = value;
                    self.readOrder.append(key);

        textFileHandle.close()
        if keep_includes:
            self.includes += fileincludes
            self.defaults += filedefaults
        else:
            for ffile in fileincludes:
                if os.path.isabs(ffile):
                    self.readFile(ffile)
                else:
                    self.readFile(os.path.join(os.path.dirname(filename), ffile))
            for ffile in filedefaults:
                if os.path.isabs(ffile):
                    self.readFile(ffile, if_not_defined=True)
                else:
                    self.readFile(os.path.join(os.path.dirname(filename), ffile), if_not_defined=True)

        return self.params

    def saveFile(self, filename):
        with open(filename, 'w') as f: f.write("\n".join(self.fileLines()))

    def fileLines(self):

        def asIniText(value):
            if type(value) == type(''): return value
            if type(value) == type(True):
                return str(value)[0]
            return str(value)

        parameterLines = []
        for include in self.includes:
            parameterLines.append('INCLUDE(' + include + ')')
        for default in self.defaults:
            parameterLines.append('DEFAULT(' + default + ')')

        keys = self.params.keys()
        keys.sort()

        for key in self.readOrder:
            if key in keys:
                parameterLines.append(key + '=' + asIniText(self.params[key]));
                keys.remove(key)
        for key in keys:
            parameterLines.append(key + '=' + asIniText(self.params[key]));

        return parameterLines


    def replaceTags(self, placeholder, text):

        for key in self.params:
            self.params[key] = self.params[key].replace(placeholder, text);

            return self.params
        
    def delete_keys(self,keys):
        for k in keys: self.params.pop(k,None)
        
    def _undefined(self,name):
        raise Exception('parameter not defined: ' + name)
        
    def bool(self,name, default=False):
        if name in self.params:
            s=self.params[name]
            if isinstance(s,bool): return s
            if s=='T': return True
            elif s=='F': return False
            raise Exception('parameter does not have valid T or F boolean value: ' + name)
        elif default is not None: return default
        else: self._undefined(name)
        
    def string(self,name, default=None):
        if name in self.params: return str(self.params[name])
        elif default is not None: return default
        else: self._undefined(name)
        
    def float(self,name,default=None):
        if name in self.params: return float(self.params[name])
        elif default is not None: return default
        else: self._undefined(name)
        
    def int(self,name,default=None):
        if name in self.params: return int(self.params[name])
        elif default is not None: return default
        else: self._undefined(name)
        
    def array_int(self,name,index=1, default=None):
        return self.int(name+'(%u)'%index,default)
        
    def array_string(self,name,index=1, default=None):
        return self.string(name+'(%u)'%index,default)

    def array_bool(self,name,index=1, default=None):
        return self.bool(name+'(%u)'%index,default)
        
    def array_float(self,name,index=1, default=None):
        return self.float(name+'(%u)'%index,default)
        