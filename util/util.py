
from datetime import datetime
import time
import sys
import os

#from callCount import FunctionCallCount

#@FunctionCallCount
def any(iterable):
    """
    any()
    Purpose:    Returns true if any element in the iterable evaluates to true otherwise, return false
    Parameters: iterable [type=list,tuple]
                    Any item that can be iterated over
    Returns:    A boolean specifying whether any elements of the iterable evaluate to true
    """
    for element in iterable:
        if element: return True
    return False


#@FunctionCallCount
def uniquify(list):
    """
    uniqify()
    Purpose:    Remove duplicate items in a list, preserving the original order
    Parameters: list [type=list]
                    A list of items (must be hashable)
    Returns:    A list of unique items
    """
    seen = {}
    result = []
    for elem in list:
        if elem in seen: continue
        seen[elem] = 1
        result.append(elem)
    return result


#@FunctionCallCount
def timeParse(time_string, time_format):
    """
    timeParse()
    Purpose:    Parse a string containing the time according to the given format
    Parameters: time_string [type=string]
                    The a string containing the time to parse.
                time_format [type=string]
                    The format to parse the time string against.
    Returns:    A datetime object set to the date/time given by time_string
    """
    time_obj = datetime(*(time.strptime(time_string, time_format)[0:6]))
    return time_obj

#@FunctionCallCount
def ravelIndex(pos, shape):
    """
    ravelIndex()
    """
    res = 0
    acc = 1
    for pi, si in zip(reversed(pos), reversed(shape)):
        res += pi * acc
        acc *= si
    return res

def warning(message):
    """
    warning()
    Purpose:    Prints a warning message.
    Parameters: message [type=string]
                    The warning message to print.
    Returns:    [nothing]
    """
    print "Warning: %s" % message
    return

def fatalError(message):
    """
    fatalError()
    Purpose:    Prints an error message and then exits the script.
    Parameters: message [type=string]
                    The error message to print.
    Returns:    [nothing]
    """
    print
    print "Error: %s" % message
    sys.exit()

def tryCreate(path):
    """
    tryCreate()
    Purpose:    Creates a directory if one does not exist already
    Parameters: path [type=string]
                    The path to Schroedinger's folder
    Returns:    [nothing]
    """
    try:
        os.stat(path)
    except:
        os.makedirs(path)

    print("Creating folder: " + path)
    return
