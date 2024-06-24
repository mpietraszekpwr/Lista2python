import itertools
from operator import mul


class Wielomian:

    def __init__(self, wsp):
        self.wsp = wsp[::-1]

    def __str__(self):
        st2 = "{}x^{}"
        st1 = "{}x"
        w = []
        for i, a in enumerate(self.wsp):
            if a != 0:
                if i > 1:
                    w.append(st2.format(a, i))
                elif i == 1:
                    w.append(st1.format(a))
                else:
                    w.append(str(a))
        return ' + '.join(w[::-1])

    def deg(self):
        return len(self.wsp) - 1

    def __add__(self, other):
        temp = itertools.zip_longest(self.wsp, other.wsp, fillvalue=0)
        n_wsp = [a + b for a, b in temp]
        return Wielomian(n_wsp[::-1])

    def __sub__(self, other):
        temp = itertools.zip_longest(self.wsp, other.wsp, fillvalue=0)
        n_wsp = [a - b for a, b in temp]
        return Wielomian(n_wsp[::-1])

    def __mul__(self, other):
        w1 = self.wsp.copy()
        w2 = other.wsp.copy()
        l = 2 * max(len(w1), len(w2))
        w1.extend([0] * (l - len(w1)))
        w2.extend([0] * (l - len(w2)))
        n_deg = self.deg() + other.deg()
        n_wsp = [sum(map(mul, w1[0:k + 1], w2[0:k + 1][::-1])) for k in range(n_deg + 1)]
        return Wielomian(n_wsp[::-1])

    def val(self, x):
        tw = self.wsp[::-1]
        v = 0
        for k in range(self.deg() + 1):
            v = (v * x + tw[k])
        return v

z1=Wielomian([2,3])
print(z1)

z1.val(2)

z2=Wielomian([2,1,-1,1])
print(z2)

z2.val(-1)

print(z1+z2)

print(z1*z2)

z1*=z2
print(z1)