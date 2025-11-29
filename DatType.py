from typing import overload,Union

class Vec_t:
    def __init__(self, x=float(0), y=float(0), z=float(0)):
        self.x = x
        self.y = y
        self.z = z
    
    @overload
    def __add__(self, other:"Vec_t"):
        ...
    
    def __add__(self, other:"Vec_t"):
        return Vec_t(self.x + other.x, self.y + other.y, self.z + other.z)
    
    @overload
    def __sub__(self, other:"Vec_t"):
        ...
    
    def __sub__(self, other:"Vec_t"):
        return Vec_t(self.x - other.x, self.y - other.y, self.z - other.z)
    
    @overload
    def __mul__(self, other:float):
        ...
    
    @overload
    def __mul__(self, other:"Vec_t"):
        ...
    
    def __mul__(self, other:Union[float, "Vec_t"]):
        if isinstance(other, float):
            return Vec_t(self.x * other, self.y * other, self.z * other)
        elif isinstance(other, Vec_t):
            return float(self.x * other.x + self.y * other.y + self.z * other.z)
        else:
            raise TypeError("Unsupported type for multiplication")

    def __truediv__(self, other:float):
        return Vec_t(self.x / other, self.y / other, self.z / other)
    
    def __str__(self):
        return f"({self.x}, {self.y}, {self.z})"
    
class Mat_t:
    def __init__(self, m00=float(0), m01=float(0), m02=float(0), m10=float(0), m11=float(0), m12=float(0), m20=float(0), m21=float(0), m22=float(0)):
        self.x0 = m00
        self.x1 = m01
        self.x2 = m02
        self.y0 = m10
        self.y1 = m11
        self.y2 = m12
        self.z0 = m20
        self.z1 = m21
        self.z2 = m22
    
    @overload
    def __mul__(self, other:float):
        ...
    
    @overload
    def __mul__(self, other:"Mat_t"):
        ...
    
    @overload
    def __mul__(self, other: Vec_t):
        ...
    
    def __mul__(self, other:Union[float, "Mat_t", Vec_t]):
        if isinstance(other, float):
            return Mat_t(self.x0 * other, self.x1 * other, self.x2 * other,
                         self.y0 * other, self.y1 * other, self.y2 * other,
                         self.z0 * other, self.z1 * other, self.z2 * other)
        elif isinstance(other, Mat_t):
            return Mat_t(self.x0 * other.x0 + self.x1 * other.y0 + self.x2 * other.z0,
                         self.x0 * other.x1 + self.x1 * other.y1 + self.x2 * other.z1,
                         self.x0 * other.x2 + self.x1 * other.y2 + self.x2 * other.z2,
                         self.y0 * other.x0 + self.y1 * other.y0 + self.y2 * other.z0,
                         self.y0 * other.x1 + self.y1 * other.y1 + self.y2 * other.z1,
                         self.y0 * other.x2 + self.y1 * other.y2 + self.y2 * other.z2,
                         self.z0 * other.x0 + self.z1 * other.y0 + self.z2 * other.z0,
                         self.z0 * other.x1 + self.z1 * other.y1 + self.z2 * other.z1,
                         self.z0 * other.x2 + self.z1 * other.y2 + self.z2 * other.z2)
        elif isinstance(other, Vec_t):
            return Vec_t(self.x0 * other.x + self.x1 * other.y + self.x2 * other.z,
                         self.y0 * other.x + self.y1 * other.y + self.y2 * other.z,
                         self.z0 * other.x + self.z1 * other.y + self.z2 * other.z)
        else:
            raise TypeError("Unsupported type for multiplication")

    def __truediv__(self, other:float):
        return Mat_t(self.x0 / other, self.x1 / other, self.x2 / other,
                     self.y0 / other, self.y1 / other, self.y2 / other,
                     self.z0 / other, self.z1 / other, self.z2 / other)

    def __add__(self, other:"Mat_t"):
        return Mat_t(self.x0 + other.x0, self.x1 + other.x1, self.x2 + other.x2,
                     self.y0 + other.y0, self.y1 + other.y1, self.y2 + other.y2,
                     self.z0 + other.z0, self.z1 + other.z1, self.z2 + other.z2)
    
    def __sub__(self, other:"Mat_t"):
        return Mat_t(self.x0 - other.x0, self.x1 - other.x1, self.x2 - other.x2,
                     self.y0 - other.y0, self.y1 - other.y1, self.y2 - other.y2,
                     self.z0 - other.z0, self.z1 - other.z1, self.z2 - other.z2)
    
    def __str__(self):
        return f"({self.x0}, {self.x1}, {self.x2} , {self.y0}, {self.y1}, {self.y2}, {self.z0}, {self.z1}, {self.z2})"

if __name__ == "__main__":
    vec = Vec_t(1, 2, 3)
    mat = Mat_t(1, 2, 3, 4, 5, 6, 7, 8, 9)

    print(vec)
    print(mat)

    print(vec * vec)
    print(vec - vec)
    print(vec + vec)

    print(mat * mat)
    print(mat + mat)
    print(mat - mat)

    print(mat * vec)