
from datetime import datetime

class FunctionCallCount(object):
    _count_hash = {}
    _time_hash = {}

    def __init__(self, func):
        self.func = func

    def __call__(self, *args, **kwargs):
        try:
            FunctionCallCount._count_hash[self.func] += 1
        except KeyError:
            FunctionCallCount._count_hash[self.func] = 1

        try:
            start = datetime.now()
            return_val = self.func(self.obj, *args, **kwargs)
            end = datetime.now()
        except AttributeError:
            start = datetime.now()
            return_val = self.func(*args, **kwargs)
            end = datetime.now()

        time_diff = end - start
        time_seconds = time_diff.seconds + float(time_diff.microseconds) / 1000000

        try:
            FunctionCallCount._time_hash[self.func] += time_seconds
        except KeyError:
            FunctionCallCount._time_hash[self.func] = time_seconds

        return return_val

    def __get__(self, instance, klass):
        self.obj = instance
        return self

def dumpCallCounts():
    for func in FunctionCallCount._count_hash.iterkeys():
        print func.__name__, FunctionCallCount._count_hash[func], FunctionCallCount._time_hash[func]

def main():
    class TestClass:
        def __init__(self):
            return

        @FunctionCallCount
        def test_member(self, arg1, arg2):
            return ", ".join([ str(i) for i in range(arg1, arg2) ])

    @FunctionCallCount
    def test_func(arg1, arg2):
        print ", ".join([ str(i) for i in range(arg1, arg2) ])

    mytest = TestClass()
    print mytest.test_member(1, 4)
    print mytest.test_member(5, 8)

    print test_func(1, 4)
    print test_func(5, 8)

    dumpCallCounts()

if __name__ == "__main__":
    main()
