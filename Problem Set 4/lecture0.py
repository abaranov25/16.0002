# 6.009 Lecture 0

import pytest

#######################################


#a = 307
#b = a
#
#print('a', a)
#print('b', b)


#######################################

x = ['baz', 302, 303, 304]
y = x

x[0] = 388

print('x', x)
print('y', y)


#######################################


#a = [301, 302, 303]
#b = [a, a, a]
#
#print('b', b)


#######################################


#x = ('baz', [301, 302], 303, 304)
#x[0] = 382
#
#print('x', x)


#######################################


#x = 500
#def foo(y):
#    return x + y
#
#z = foo(307)
#
#print('z', z)


#######################################


#def bar(x):
#    x = 2000
#    return foo(307)
#
#z = bar(1000)
#
#print('z', z)


#######################################

functions = []
for i in range(5):
    def func(x):
        return x + i
    functions.append(func)

for f in functions:
    print(f(12))