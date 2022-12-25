from constraints.constraint import Constraint


class Horizontal(Constraint):
    """Класс горизонтальности"""
    def __init__(self, Points):
        """Constructor"""
        self.Points = Points
        self.Ls = 1
        self.Name = "Горизонтальность"

    def LocalCon(self, D, L):
        a = 5
        matrix = [0] * a
        for i in range(a):
            matrix[i] = [0] * a


        matrix[0][0] = 1
        matrix[1][1] = 1
        matrix[2][2] = 1
        matrix[3][3] = 1

        matrix[1][4] = -1
        matrix[3][4] = 1
        matrix[4][1] = -1
        matrix[4][3] = 1



        F = [0] * 5

        F[0] = -(D[0][0])
        F[1] = -(D[0][1] - L[0])
        F[2] = -(D[1][0])
        F[3] = -(D[1][1] + L[0])

        F[4] = -(self.Points[1].aCoord[1]+D[1][1] - self.Points[0].aCoord[1] - D[0][1])

        return matrix, F
