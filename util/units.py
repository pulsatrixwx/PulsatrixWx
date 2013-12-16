
import re
from copy import copy

from util import any, fatalError
#from callCount import FunctionCallCount

class Units(object):
    """
    Units
    Purpose:    Handles the parsing of unit strings, evaluation of units functions, and conversions of values from one units system to another.
    Started:    2 February 2011 by Tim Supinie (tsupinie@ou.edu)
    Completed   [not yet]
    Modified:   [not yet]
    """
    # List of SI units
    _metric_units = ['m', 's', 'g', 'K', 'DK', 'J', 'Pa','W']

    # Dictionary of SI prefixes to their values.  Conversions are done by multiplying by (value_for_from_prefix / value_for_to_prefix) ** exponent
    _metric_prefixes = {
        'p':1e-12,
        'n':1e-9,
        'u':1e-6,
        'm':1e-3,
        'c':1e-2,
        'd':1e-1,
        '':1e0,
        'da':1e1,
        'h':1e2,
        'k':1e3,
        'M':1e6,
        'G':1e9,
        'T':1e12,
    }

    # List of units in the imperial system (okay, really, they're just units that aren't strictly SI units).
    _imperial_units = ['F', 'DF', 'C', 'DC', 'inHg', 'in', 'b', 'NM', 'mi', 'ft', 'min', 'hr', 'dy', '%', '', 'dBZ']

    # Conversion functions between imperial and metric units systems.  The structure of is a dictionary of dictionaries.  The top-level dictionary is
    #   the units "from" and the second-level dictionaries are the units "to."
    _imperial_conversions = {
        # Temperature 
       'C':{
            'K':lambda(temp_C): temp_C + 273.15,
            'F':lambda(temp_C): 1.8 * temp_C + 32,
        },
        'K':{
            'C':lambda(temp_K): temp_K - 273.15,
            'F':lambda(temp_K): 1.8 * temp_K - 459.67,
        },
        'F':{
            'C':lambda(temp_F): 5 * (temp_F - 32) / 9.,
            'K':lambda(temp_F): 5 * (temp_F + 459.67) / 9.,
        },

        # Temperature Change
        'DC':{
            'DK':lambda(temp_change_C): temp_change_C,
            'DF':lambda(temp_change_C): temp_change_C * 1.8,
        },
        'DK':{
            'DC':lambda(temp_change_K): temp_change_K,
            'DF':lambda(temp_change_K): temp_change_K * 1.8,
        },
        'DF':{
            'DC':lambda(temp_change_F): temp_change_F / 1.8,
            'DK':lambda(temp_change_F): temp_change_F / 1.8,
        },

        # Pressure
        'inHg':{
            'b':lambda(pres_inHg): pres_inHg * .033863958,
            'Pa':lambda(pres_inHg): pres_inHg * 3386.3958,
        },
        'b':{
            'Pa':lambda(pres_b): pres_b * 1e5,
            'inHg':lambda(pres_b): pres_b * 29.5250192,
        },
        'Pa':{
            'b':lambda(pres_Pa): pres_Pa * 1e-5,
            'inHg':lambda(pres_Pa): pres_Pa * 2.9520192e-4,
        },

        #Length
        'm':{
            'in':lambda(len_m): len_m * 39.3700787,
            'ft':lambda(len_m): len_m * 3.2808399,
            'mi':lambda(len_m): len_m / 1609.344,
            'NM':lambda(len_m): len_m / 1852.,
        },
        'ft':{
            'm':lambda(len_ft): len_ft / 3.2808399,
            'mi':lambda(len_ft): len_ft / 5280.,
            'NM':lambda(len_ft): len_ft / 6067.11549,
            'in':lambda(len_ft): len_ft * 12.0,
        },
        'mi':{
            'm':lambda(len_mi): len_mi * 1609.344,
            'ft':lambda(len_mi): len_mi * 5280.,
            'NM':lambda(len_mi): len_mi / 1.15077945,
        },
        'NM':{
            'm':lambda(len_NM): len_NM * 1852,
            'mi':lambda(len_NM): len_NM * 1.15077945,
            'ft':lambda(len_NM): len_NM * 6076.11549,
        },

    #Time
        's':{
            'min':lambda(time_s): time_s / 60.,
            'hr':lambda(time_s): time_s / 3600.,
            'dy':lambda(time_s): time_s / (24. * 3600.),
        },
        'min':{
            's':lambda(time_min): time_min * 60,
            'hr':lambda(time_min): time_min / 60.,
            'dy':lambda(time_min): time_min / (24. * 60.),
        },
        'hr':{
            's':lambda(time_hr): time_hr * 3600,
            'min':lambda(time_hr): time_hr * 60,
            'dy': lambda(time_hr): time_hr / 24.,
        },
        'dy':{
            'hr':lambda(time_dy): time_dy * 24,
            'min':lambda(time_dy): time_dy * 24 * 60,
            's':lambda(time_dy): time_dy * 24 * 3600,
        },

    # Fractional
        '%': {
            '':lambda(frac_percent): frac_percent / 100.
        },
        '': {
            '%':lambda(frac_unitless): frac_unitless * 100.
        }
    }

    # Unit synonyms
    _replacements = {
        'kts':'NM hr-1'
    }

    _replacements_inv = dict(zip(_replacements.values(), _replacements.keys()))    

    # List of operators
    _operators = ["+", "-", "*", "/"]

    def __init__(self, string):
        """
        __init__()
        Purpose:    Constuctor for the Units class (parses unit strings and stores them internally as a list of dictionaries).
        Parameters: string [type=str]
                        Units string to parse and store.
        """
        self._scalar = 1
        bit_re = re.compile("^(%s)?(%s)([\\d\\-]+)?$" % ("|".join(Units._metric_prefixes), "|".join(["|".join(Units._metric_units), "|".join(Units._imperial_units)])))
        unit_bits = re.split("[\\s]+", string)
        self._units = []


        if string == "":
            self._units.append({'pre':"", 'unt':"", 'exp':1})
        else:
            for bit in unit_bits:
                if bit in Units._replacements:
                    replacement = Units(Units._replacements[bit])
                    self._units.extend(replacement._units)
                else:
                    try:
                        self._scalar = float(bit)
                        continue
                    except ValueError:
                        try:
                            self._scalar = int(bit)
                            continue
                        except ValueError:
                            pass

                    match = bit_re.search(bit)
                    try:
                        prefix, unit, exponent = match.groups(default="")
                    except AttributeError:
                        fatalError("Bit '%s' in unit string '%s' could not be parsed." % (bit, string))
            
                    self._units.append({'pre':prefix, 'unt':unit, 'exp':int(exponent) if exponent != "" else 1})
        return

    def __str__(self):
        """
        __str__() [public]
        Purpose:    Returns a string representation of the Units object (as used by 'print', str(), and others).
        Parameters: [nothing]
        Returns:    A units string containing the same units as the object
        """
        strings = []
        for u in self._units:
            pre, unt, exp = u.values()

            if exp == 1:
                exp = ""
            else:
                exp = str(exp)
            strings.append("%s%s%s" % (pre, unt, exp))

        if self._scalar != 1:
            return "%1.e %s" % (self._scalar, " ".join(strings))
        else:
            return " ".join(strings)

    def toLaTexString(self):
        """
        toLatexString() [public]
        Purpose:    Returns a string representation of the unit with special characters in LaTex.
        Parameters: [nothing]
        Returns:    A string representation of the units with LaTex special characters.
        """
        strings = []

        if str(self) in Units._replacements_inv:
            return Units._replacements_inv[str(self)]

        for u in self._units:
            pre, unt, exp = u.values()

            if exp == 1:
                exp = ""
            else:
                exp = "$^{%s}$" % exp

            if pre == "u":
                pre = r"$\mathrm{\mu}$"

            if unt in [ 'F', 'C' ]:
                #unt = r"$^\circ$%s" % unt
                unt = unt
            elif unt in [ 'DF', 'DC' ]:
                unt = r"$\mathrm{\Delta}^{\circ}$%s" % unt[1:] 
            elif unt == 'DK':
                unt = r"$\mathrm{\Delta}K"

            strings.append("%s%s%s" % (pre, unt, exp))

        if self._scalar != 1:
            scalar_exp = ("%1.e" % self._scalar).split('e+')
            if scalar_exp[0] == "1":
                return "$10^{%s}$ %s" % (scalar_exp[1], " ".join(strings))
            else:
                return "$%s \times 10^{%s}$ %s" % (scalar_exp[0],scalar_exp[1], " ".join(strings))
        else:
            return " ".join(strings)

    def __contains__(self, other):
        """
        __contains__() [public]
        Purpose:    Check to see whether this Units object contains a given object (as used by 'in').
        Parameters: other [type=str,dict,Units]
                        The units object to check for.  It can be either a string (in which case it will be parsed to a Units object), dictionary
                        containing the internal representation of the units object, or a Units object.
        Returns:    A boolean specifying whether or not the other object is contained in this object.
        """
        if type(other) == str:
            other = Units(other)

        if len(other) > 1:
            pass

        is_equal = False

        if type(other) == dict:
            other_pre = other['pre']
            other_unt = other['unt']
            other_exp = other['exp']
        else:
            other_pre = other._units[0]['pre']
            other_unt = other._units[0]['unt']
            other_exp = other._units[0]['exp']

        for bit in self._units:
            is_equal |= (bit['pre'] == other_pre and bit['unt'] == other_unt and bit['exp'] == other_exp)
        return is_equal

    def __eq__(self, other):
        """
        __eq__() [public]
        Purpose:    Check to see if another string is equal to this string (contains the exact same units).
        Parameters: other [type=str,Units]
                        The units object to test.  It can be either a string (in which case it will be parsed to a Units object), or a Units object.
        Returns:    A boolean specifying whether or not the two units are equal.
        """
        if type(other) == str:
            other = Units(other)

        is_equal = other._scalar == self._scalar

        for bit in self._units:
            is_equal &= bit in other
        return is_equal

    def __ne__(self, other):
        return not self.__eq__(other)

    def __len__(self):
        return len(self._units)

    def __add__(self, other):
        if type(other) == str:
            other = Units(other)

        if self != other:
            fatalError("Can't add units '%s' and '%s'." % (self, other))
        return copy(self)

    def __iadd__(self, other):
        if type(other) == str:
            other = Units(other)

        if self != other:
            fatalError("Can't add units '%s' and '%s'." % (self, other))
        return self

    def __sub__(self, other):
        if type(other) == str:
            other = Units(other)

        if self != other and not (len(other) == 1 and other._units[0]['unt'] in ['DK', 'DC', 'DF']):
            fatalError("Can't subtract units '%s' and '%s'." % (self, other))
        result = copy(self)

        if len(result) == 1 and result._units[0]['unt'] in ['K', 'C', 'F']:
            result._units[0]['unt'] = "D%s" % result._units[0]['unt']
        elif len(other) == 1 and other._units[0]['unt'] in ['DK', 'DC', 'DF']:
            result._units[0]['unt'] = result._units[0]['unt'][-1]

        return result

    def __isub__(self, other):
        if type(other) == str:
            other = Units(other)

        if self != other and not (len(other) == 1 and other._units[0]['unt'] in ['DK', 'DC', 'DF']):
            fatalError("Can't subtract units '%s' and '%s'." % (self, other))

        if len(other) == 1 and other._units[0]['unt'] in ['K', 'C', 'F']:
            self._units[0]['unt'] = "D%s" % self._units[0]['unt']
        elif len(other) == 1 and other._units[0]['unt'] in ['DK', 'DC', 'DF']:
            self._units[0]['unt'] = self._units[0]['unt'][-1]

        return self

    def __mul__(self, other):
        result = copy(self)
        result._units.extend(other._units)
        result._scalar *= other._scalar
        result._simplify()
        return result

    def __imul__(self, other):
        self._units.extend(other._units)
        self._scalar *= other._scalar
        self._simplify()
        return self

    def __div__(self, other):
        result = copy(self)
        for bit in other._units:
            result._units.append({'pre':bit['pre'], 'unt':bit['unt'], 'exp':-bit['exp']})

        result._scalar /= other._scalar

        result._simplfy()
        return result

    def __idiv__(self, other):
        for bit in other._units:
            self._units.append({'pre':bit['pre'], 'unt':bit['unt'], 'exp':-bit['exp']})

        self._scalar /= other._scalar

        self._simplify()
        return self

    def __getitem__(self, index):
        return self._units[index]

    def _simplify(self):
        simplifications = []
        for idx in xrange(len(self._units) - 1):
            for jdy in xrange(idx + 1, len(self._units)):
                if self._units[idx]['pre'] == self._units[jdy]['pre'] and self._units[idx]['unt'] == self._units[jdy]['unt']:
                    simplifications.append((idx, jdy))

        for idx, jdy in simplifications:
            self._units[idx]['exp'] += self._units[jdy]['exp']

        if len(simplifications) > 0:
            idxs, jdys = zip(*simplifications)
            jdys = sorted(jdys, reverse=True)

            for jdy in jdys:
                del self._units[jdy]

        idx = 0
        while idx < len(self._units):
            if self._units[idx]['exp'] == 0:
                del self._units[idx]
            else:
                idx += 1

        return

    #@FunctionCallCount
    @classmethod
    def convert(klass, value, units_from, units_to, force_exception=False):
        """
        convert() [public, static]
        Purpose:    Convert a value from one unit to another
        Parameters: value [type=int,float,np.array]
                        The value or array of values to convert
                    units_from [type=string]
                        The units of the incoming value(s)
                    units_to [type=string]
                        The units to which to convert the value(s)
        Returns:    The value in the new unit system.
        """
        if type(units_from) == str:
            units_from = Units(units_from)

        if type(units_to) == str:
            units_to = Units(units_to)

        if len(units_from) != len(units_to):
            fatalError("Cannot convert between units '%s' and '%s' (units are of different lengths)." % (units_from, units_to))

        for idx in xrange(len(units_from)):
            u_from = units_from[idx]
            u_to = units_to[idx]

            if u_from['unt'] != u_to['unt']:
                # Only attempt converting the base unit if they're different
                for idx in xrange(abs(u_from['exp'])):
                    # Do the conversion as many times as the exponent (e.g. m^2 to ft^2, need to multiply by 3.281 twice)
                    try:
                        # If the exponent is greater than 0, convert from the old to the new, otherwise, do the inverse.
                        if u_from['exp'] > 0: value = Units._imperial_conversions[u_from['unt']][u_to['unt']](value)
                        else:                 value = Units._imperial_conversions[u_to['unt']][u_from['unt']](value)
                    except KeyError:
                        # The from and to units are not in the same section of the nested dictionary, so we're trying to do something like convert
                        #   m to s or something like that.  Error.
                        if not force_exception:
                            fatalError("Cannot convert units from '%s' to '%s'." % (units_from, units_to))
                        else:
                            raise ValueError("Cannot convert units from '%s' to '%s'" % (units_from, units_to))
            # Now that we know the value is in the same unit system, do any metric prefix conversions.
            value = value * (Units._metric_prefixes[u_from['pre']] / Units._metric_prefixes[u_to['pre']]) ** u_from['exp']

        return value

    @classmethod
    def quantity(klass, units):

        if units == "%":
            # Check for units of percent
            return "percent"
        elif units == "":
            # Check for no units
            return "unitless"

        try:
            # Check for length units
            Units.convert(0, units, 'm', force_exception=True)
            return "length"
        except ValueError:
            pass

        try:
            # Check for absolute temperature units 
            Units.convert(0, units, 'K', force_exception=True)
            return "absolute-temperature"
        except ValueError:
            pass

        try:
            # Check for relative temperature units
            Units.convert(0, units, 'DK', force_exception=True)
            return "relative-temperature"
        except ValueError:
            pass

        try:
            # Check for time units
            Units.convert(0, units, 's', force_exception=True)
            return "time"
        except ValueError:
            pass

        try:
            # Check for mass units
            Units.convert(0, units, 'g', force_exception=True)
            return "mass"
        except ValueError:
            pass

        try:
            # Check for pressure units
            Units.convert(0, units, 'Pa', force_exception=True)
            return "pressure"
        except ValueError:
            pass

        # Eh? :/
        return "unknown"


    @classmethod
    def _evaluateSubFunction(klass, subfunction):
        operator_re = re.compile(r"(?<=[\D])(?:%s)(?=[\D])" % "|".join([ r"\%s" % o for o in Units._operators ]))
        unit_strings = [ Units(s.strip()) for s in operator_re.split(subfunction) ]
        operators = operator_re.findall(subfunction)

        for op in reversed(Units._operators):
            while True:
                try:
                    op_index = operators.index(op)
                except ValueError:
                    break

                exec "unit_strings[op_index] %s= unit_strings[op_index + 1]" % op 
                del operators[op_index]
                del unit_strings[op_index + 1]

        return str(unit_strings[0])

    #@FunctionCallCount
    @classmethod
    def evaluateFunction(klass, function):
        """
        evaluateFunction()
        Purpose:    Take a combination of unit strings joined by mathematical operators (a "unit function") and, if possible, collapse it into a single unit string.
        Parameters: function [type=str]
                        The unit function string to evaluate.
        Returns:    The unit string that results from the evaluation of the function.
        """
        parenthesis_re = re.compile(r"\(([^\(\)]+)\)")

        wrap = False
        if function.find(",") > -1:
            wrap = True

        while True:
            sub_functions = parenthesis_re.findall(function)

            if sub_functions == []:
                break

            for sub_func in sub_functions:
                vec_units = re.split(r"\,[\s]+", sub_func)
                eval_sub_func = sub_func
                for vec in vec_units:
                    eval_vec = Units._evaluateSubFunction(vec)
                    eval_sub_func = eval_sub_func.replace(vec, eval_vec)
                function = function.replace("(%s)" % sub_func, eval_sub_func)

        if wrap:
            func_parts = re.split(r"\,[\s]+", function)
            all_equal = True
            for idx in xrange(len(func_parts) - 1):
                for jdy in xrange(idx + 1, len(func_parts)):
                    all_equal &= func_parts[idx] == func_parts[jdy]

            if not all_equal:
                fatalError("Vector components do not have equal units.")

            function = func_parts[0]
        else:
            function = Units._evaluateSubFunction(function)

        return function

if __name__ == "__main__":
    u1 = Units("kts")
    u2 = Units("K")
    print u1
    print u2
    print u1 == u2

    print Units.convert(300, "K", "C"), "(should be 26.85)"
    print Units.convert(50, "F", "C"), "(should be 10.0)"
    print Units.convert(15, "m s-1", "kts"), "(should be 29.1576674)"
    print Units.convert(15, "m s-1", "km hr-1"), "(should be 54.0)"
    print Units.convert(1500, "m", "dam"), "(should be 150.0)"

    print Units.evaluateFunction("(m2 s-2 + m2 s-2) / (m s-2 + m s-2) - m")
    print Units.evaluateFunction("(m s-1, m s-1)")
    print Units.evaluateFunction("(m / s, m / s)")
