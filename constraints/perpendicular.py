from constraints.constraint import Constraint


class Perpendicular(Constraint):
    """Класс перпендикулярности"""
    def __init__(self, Points):
        """Constructor"""
        self.Points = Points
        self.Ls = 1
        self.Name = "Перпендикулярность"

    def get_description(self):
        return "Перпендикулярность отрезков {" + str(self.Points[0].v_return()) + "; " + str(self.Points[1].v_return()) \
               + "} и {" + str(self.Points[2].v_return()) + "; " + str(self.Points[3].v_return()) + "}"

    def LocalCon(self, deltas, lambdas):
        k = 9
        matrix = [0] * k
        for i in range(k):
            matrix[i] = [0] * k

        a = self.Points[1].aCoord[0] + deltas[1][0] - self.Points[0].aCoord[0] - deltas[0][0]
        b = self.Points[1].aCoord[1] + deltas[1][1] - self.Points[0].aCoord[1] - deltas[0][1]
        c = self.Points[3].aCoord[0] + deltas[3][0] - self.Points[2].aCoord[0] - deltas[2][0]
        d = self.Points[3].aCoord[1] + deltas[3][1] - self.Points[2].aCoord[1] - deltas[2][1]

        lam = lambdas[0]

        for i in range(8):
            matrix[i][i] = 1

        matrix[8][0] = matrix[0][8] = -c
        matrix[8][1] = matrix[1][8] =-d
        matrix[8][2] = matrix[2][8] = c
        matrix[8][3] = matrix[3][8] = d
        matrix[8][4] = matrix[4][8] = -a
        matrix[8][5] = matrix[5][8] = -b
        matrix[8][6] = matrix[6][8] = a
        matrix[8][7] = matrix[7][8] = b

        F = [0] * k

        F[0] = -(deltas[0][0] - lam * c)
        F[1] = -(deltas[0][1] - lam * d)
        F[2] = -(deltas[1][0] + lam * c)
        F[3] = -(deltas[1][1] + lam * d)
        F[4] = -(deltas[2][0] - lam * a)
        F[5] = -(deltas[2][1] - lam * b)
        F[6] = -(deltas[3][0] + lam * a)
        F[7] = -(deltas[3][1] + lam * b)

        F[8] = -(a*c + b*d)

        return matrix, F
