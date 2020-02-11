import os
import numpy as np


class IniError(Exception):
    pass


class IniFile:
    """
    Class for storing option parameter values and reading/saving to file

    Unlike standard .ini files, IniFile allows inheritance, in that a .ini file can use
    INCLUDE(..) and DEFAULT(...) to include or override settings in another file (to avoid duplication)

    :ivar params: dictionary of name, values stored
    :ivar comments: dictionary of optional comments for parameter names
    """

    def __init__(self, settings=None, keep_includes=False, expand_environment_variables=True):
        """

        :param settings: a filename of a .ini file to read, or a dictionary of name/values
        :param keep_includes:
             - False: load all INCLUDE and DEFAULT files, making one params dictionary
             - True: only load settings in main file, and store INCLUDE and DEFAULT entries into defaults
               and includes filename lists.
        :param expand_environment_variables: whether to expand $(var) placeholders in parameter values
               using environment variables
        """

        self.params = dict()
        self.comments = dict()
        self.readOrder = []
        self.defaults = []
        self.includes = []
        self.original_filename = None
        self.expand_environment_variables = expand_environment_variables
        if isinstance(settings, str):
            self.readFile(settings, keep_includes)
        elif settings:
            self.params.update(settings)

    def expand_placeholders(self, s):
        """Expand shell variables of the forms $(var), like in Makefiles"""
        if '$(' not in s:
            return s
        res = ''
        index = 0
        pathlen = len(s)
        while index < pathlen:
            c = s[index]
            if c == '$':
                if s[index + 1] == '$':
                    res = res + c
                    index += 1
                elif s[index + 1] == '(':
                    s = s[index + 2:]
                    pathlen = len(s)
                    index = s.index(')')
                    var = s[:index]
                    if var in os.environ:
                        res = res + os.environ[var]
            else:
                res = res + c
            index += 1
        return res

    def readFile(self, filename, keep_includes=False, if_not_defined=False):
        try:
            fileincludes = []
            filedefaults = []
            self.original_filename = filename
            comments = []
            with open(filename, encoding='utf-8-sig') as textFileHandle:
                # Remove blank lines and comment lines from the python list of lists.
                for line in textFileHandle:
                    s = line.strip()
                    if s == 'END':
                        break
                    if s.startswith('#'):
                        comments.append(s[1:].rstrip())
                        continue
                    elif s.startswith('INCLUDE('):
                        fileincludes.append(s[s.find('(') + 1:s.rfind(')')])
                    elif s.startswith('DEFAULT('):
                        filedefaults.append(s[s.find('(') + 1:s.rfind(')')])
                    elif s != '':
                        eq = s.find('=')
                        if eq >= 0:
                            key = s[0:eq].strip()
                            if key in self.params:
                                if if_not_defined:
                                    continue
                                raise IniError('Error: duplicate key: ' + key + ' in ' + filename)
                            value = s[eq + 1:].strip()
                            if self.expand_environment_variables:
                                value = self.expand_placeholders(value)
                            self.params[key] = value
                            self.readOrder.append(key)
                            if len(comments):
                                self.comments[key] = comments
                    if not s.startswith('#'):
                        comments = []

            if keep_includes:
                self.includes += fileincludes
                self.defaults += filedefaults
            else:
                for ffile in fileincludes:
                    if os.path.isabs(ffile):
                        self.readFile(ffile, if_not_defined=if_not_defined)
                    else:
                        self.readFile(os.path.join(os.path.dirname(filename), ffile), if_not_defined=if_not_defined)
                for ffile in filedefaults:
                    if os.path.isabs(ffile):
                        self.readFile(ffile, if_not_defined=True)
                    else:
                        self.readFile(os.path.join(os.path.dirname(filename), ffile), if_not_defined=True)

            return self.params
        except:
            print('Error in ' + filename)
            raise

    def __str__(self):
        return "\n".join(self.fileLines())

    def saveFile(self, filename=None):
        """
        Save to a .ini file

        :param filename:  name of file to save to
        """

        if not filename:
            filename = self.original_filename
        if not filename:
            raise IniError('No filename for iniFile.saveFile()')
        with open(filename, 'w', encoding='utf-8') as f:
            f.write(str(self))

    def fileLines(self):

        def asIniText(value):
            if isinstance(value, type('')):
                return value
            if type(value) == bool:
                return str(value)[0]
            return str(value)

        parameterLines = []
        for include in self.includes:
            parameterLines.append('INCLUDE(' + include + ')')
        for default in self.defaults:
            parameterLines.append('DEFAULT(' + default + ')')

        keys = list(self.params.keys())
        keys.sort()

        for key in self.readOrder:
            if key in keys:
                parameterLines.append(key + '=' + asIniText(self.params[key]))
                keys.remove(key)
        for key in keys:
            parameterLines.append(key + '=' + asIniText(self.params[key]))

        return parameterLines

    def replaceTags(self, placeholder, text):
        for key in self.params:
            self.params[key] = self.params[key].replace(placeholder, text)
        return self.params

    def delete_keys(self, keys):
        for k in keys:
            self.params.pop(k, None)

    def _undefined(self, name):
        raise IniError('parameter not defined: ' + name)

    def hasKey(self, name):
        """
        Test if key name exists

        :param name: parameter name
        :return: True or False test if key name exists
        """
        return name in self.params

    def isSet(self, name, allowEmpty=False):
        """
        Tests whether value for name is set or is empty

        :param name: name of parameter
        :param allowEmpty: whether to allow empty strings (return True is parameter name exists but is not set, "x = ")
        """

        return name in self.params and (allowEmpty or self.params[name] != "")

    def asType(self, name, tp, default=None, allowEmpty=False):
        if self.isSet(name, allowEmpty):
            if tp == bool:
                return self.bool(name, default)
            elif tp == list:
                return self.split(name, default)
            elif tp == np.ndarray:
                return self.ndarray(name, default)
            else:
                return tp(self.params[name])
        elif default is not None:
            return default
        else:
            self._undefined(name)

    def setAttr(self, name, instance, default=None, allowEmpty=False):
        """
        Set attribute of an object to value of parameter, using same type as existing value or default

        :param name: parameter name
        :param instance: instance of an object, so instance.name is the value to set
        :param default: default value if instance.name does not exist
        :param allowEmpty: whether to allow empty values
        """
        default = getattr(instance, name, default)
        setattr(instance, name, self.asType(name, type(default), default, allowEmpty=allowEmpty))

    def getAttr(self, instance, name, default=None, comment=None):
        val = getattr(instance, name, default)
        self.params[name] = val
        if comment:
            self.comments[name] = comment

    def bool(self, name, default=False):
        """
        Get boolean value

        :param name: parameter name
        :param default: default value if not set
        """
        if self.isSet(name):
            s = self.params[name]
            if isinstance(s, bool):
                return s
            if s[0] == 'T':
                return True
            elif s[0] == 'F':
                return False
            raise IniError('parameter does not have valid T(rue) or F(alse) boolean value: ' + name)
        elif default is not None:
            return default
        else:
            self._undefined(name)

    def bool_list(self, name, default=None):
        """
        Get list of boolean values, e.g. from name = T F T

        :param name: parameter name
        :param default: default value if not set
        """

        if not default:
            default = []
        return self.split(name, default, tp=bool)

    def string(self, name, default=None, allowEmpty=True):
        """
        Get string value

        :param name: parameter name
        :param default: default value if not set
        :param allowEmpty: whether to return empty string if value is empty (otherwise return default)
        """
        return self.asType(name, str, default, allowEmpty=allowEmpty)

    def list(self, name, default=None, tp=None):
        """
        Get list (from space-separated values)

        :param name: parameter name
        :param default: default value
        :param tp: type for each member of the list
        """

        if not default:
            default = []
        return self.split(name, default, tp)

    def float(self, name, default=None):
        """
        Get float value

        :param name: parameter name
        :param default: default value
        """
        return self.asType(name, float, default)

    def float_list(self, name, default=None):
        """
        Get list of float values

        :param name: parameter name
        :param default: default value
        """

        if not default:
            default = []
        return self.split(name, default, tp=float)

    def int(self, name, default=None):
        """
        Get int value

        :param name: parameter name
        :param default: default value
        """

        return self.asType(name, int, default)

    def int_list(self, name, default=None):
        """
        Get list of int values

        :param name: parameter name
        :param default: default value
        """

        if not default:
            default = []
        return self.split(name, default, tp=int)

    def split(self, name, default=None, tp=None):
        """
        Gets a list of values, optionally cast to type tp

        :param name: parameter name
        :param default: default value
        :param tp: type for each list member
        """
        if name in self.params and isinstance(self.params[name], (list, tuple)):
            if tp is None:
                return self.params[name]
            else:
                return [tp(x) for x in self.params[name]]

        s = self.string(name, default)
        if isinstance(s, str):
            if tp is not None:
                return [tp(x) for x in s.split()]
            return s.split()
        else:
            return s

    def ndarray(self, name, default=None, tp=np.float64):
        """
        Get numpy array of values

        :param name: parameter name
        :param default: default value
        :param tp: type for array
        """
        return np.array(self.split(name, default, tp=tp))

    def array_int(self, name, index=1, default=None):
        """
        Get one int value, for entries of the form name(index)

        :param name: base parameter name
        :param index: index (in brackets)
        :param default: default value
        """
        return self.int(name + '(%u)' % index, default)

    def array_string(self, name, index=1, default=None):
        """
        Get one str value, for entries of the form name(index)

        :param name: base parameter name
        :param index: index (in brackets)
        :param default: default value
        """

        return self.string(name + '(%u)' % index, default)

    def array_bool(self, name, index=1, default=None):
        """
        Get one boolean value, for entries of the form name(index)

        :param name: base parameter name
        :param index: index (in brackets)
        :param default: default value
        """

        return self.bool(name + '(%u)' % index, default)

    def array_float(self, name, index=1, default=None):
        """
        Get one float value, for entries of the form name(index)

        :param name: base parameter name
        :param index: index (in brackets)
        :param default: default value
        """

        return self.float(name + '(%u)' % index, default)

    def relativeFileName(self, name, default=None):
        s = self.string(name, default)
        if not os.path.isabs(s) and self.original_filename is not None:
            return os.path.join(os.path.dirname(self.original_filename), s)
        return s
