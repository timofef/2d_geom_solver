import math

from constraints.constraint import Constraint


class Angle(Constraint):
    """Класс угла"""
    def __init__(self, Points, angle):
        """Constructor"""
        self.Points = Points
        self.Ls = 1
        self.Name = "Угол"
        self.u = math.radians(angle)

    def get_description(self):
        return "Угол между отрезками {" + str(self.Points[0].v_return()) + "; " + str(self.Points[1].v_return()) \
               + "} и {" + str(self.Points[2].v_return()) + "; " + str(self.Points[3].v_return()) + "}"

    def LocalCon(self, Del, lambdas):
        k = 9
        matrix = [0] * k
        for i in range(k):
            matrix[i] = [0] * k

        a = self.Points[1].aCoord[0] + Del[1][0] - self.Points[0].aCoord[0] - Del[0][0]
        b = self.Points[1].aCoord[1] + Del[1][1] - self.Points[0].aCoord[1] - Del[0][1]

        c = self.Points[3].aCoord[0] + Del[3][0] - self.Points[2].aCoord[0] - Del[2][0]
        d = self.Points[3].aCoord[1] + Del[3][1] - self.Points[2].aCoord[1] - Del[2][1]

        lam = lambdas[0]

        A = (a*c*c + b*c*d) - a*(c*c + d*d)*math.cos(self.u)*math.cos(self.u)
        B = (b*d*d + a*c*d) - b*(c*c + d*d)*math.cos(self.u)*math.cos(self.u)

        C = (c*a*a + a*b*d) - c*(a*a + b*b)*math.cos(self.u)*math.cos(self.u)
        D = (d*b*b + a*b*c) - d*(a*a + b*b)*math.cos(self.u)*math.cos(self.u)

        for i in range(8):
            matrix[i][i] = 1

        matrix[0][8] = matrix[8][0] = -2 * A
        matrix[1][8] = matrix[8][1] = -2 * B
        matrix[2][8] = matrix[8][2] = 2 * A
        matrix[3][8] = matrix[8][3] = 2 * B
        matrix[4][8] = matrix[8][4] = -2 * C
        matrix[5][8] = matrix[8][5] = -2 * D
        matrix[6][8] = matrix[8][6] = 2 * C
        matrix[7][8] = matrix[8][7] = 2 * D

        F = [0] * k

        F[0] = -(Del[0][0] + 2*lam*(-A))
        F[1] = -(Del[0][1] + 2*lam*(-B))

        F[2] = -(Del[1][0] + 2*lam*(A))
        F[3] = -(Del[1][1] + 2*lam*(B))

        F[4] = -(Del[2][0] + 2*lam*(-C))
        F[5] = -(Del[2][1] + 2*lam*(-D))

        F[6] = -(Del[3][0] + 2*lam*(C))
        F[7] = -(Del[3][1] + 2*lam*(D))

        F[8] = -(a*a*c*c + 2*a*b*c*d + b*b*d*d - (a*a+b*b)*(c*c+d*d)*math.cos(self.u)*math.cos(self.u))

        return matrix, F
