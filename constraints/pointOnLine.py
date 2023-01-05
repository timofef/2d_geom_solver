from constraints.constraint import Constraint


class PointOnLine(Constraint):
    """Класс точки на прямой"""
    def __init__(self, Points):
        """Constructor"""
        self.Points = Points
        self.Ls = 1
        self.Name = "Точка на прямой"

    def get_description(self):
        return "Точка " + str(self.Points[2].v_return()) + " на прямой {" + str(self.Points[0].v_return()) + "; " + str(self.Points[1].v_return()) + "}"

    def LocalCon(self, deltas, lambdas):
        k = 7
        matrix = [0] * k
        for i in range(k):
            matrix[i] = [0] * k

        a = self.Points[2].aCoord[0] + deltas[2][0] - self.Points[0].aCoord[0] - deltas[0][0]
        b = self.Points[2].aCoord[1] + deltas[2][1] - self.Points[0].aCoord[1] - deltas[0][1]
        c = self.Points[1].aCoord[0] + deltas[1][0] - self.Points[2].aCoord[0] - deltas[2][0]
        d = self.Points[1].aCoord[1] + deltas[1][1] - self.Points[2].aCoord[1] - deltas[2][1]

        lam = lambdas[0]

        for i in range(6):
            matrix[i][i] = 1

        matrix[6][0] = matrix[0][6] = -d
        matrix[6][1] = matrix[1][6] = c
        matrix[6][2] = matrix[2][6] = -b
        matrix[6][3] = matrix[3][6] = a
        matrix[6][4] = matrix[4][6] = b + d
        matrix[6][5] = matrix[5][6] = -a - c

        F = [0] * k

        F[0] = -(deltas[0][0] - lam * d)
        F[1] = -(deltas[0][1] + lam * c)

        F[2] = -(deltas[1][0] - lam * b)
        F[3] = -(deltas[1][1] + lam * a)

        F[4] = -(deltas[2][0] + lam * (d + b))
        F[5] = -(deltas[2][1] - lam * (c + a))

        F[6] = -(a*d - c*b)
        return matrix, F
