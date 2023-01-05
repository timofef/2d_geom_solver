from constraints.constraint import Constraint


class Horizontal(Constraint):
    """Класс горизонтальности"""
    def __init__(self, Points):
        """Constructor"""
        self.Points = Points
        self.Ls = 1
        self.Name = "Горизонтальность"

    def get_description(self):
        return "Горизонтальность отрезка {" + str(self.Points[0].v_return()) + "; " + str(self.Points[1].v_return()) + "}"

    def LocalCon(self, D, L):
        a = 5
        matrix = [0] * a
        for i in range(a):
            matrix[i] = [0] * a


        for i in range(4):
            matrix[i][i] = 1

        matrix[1][4] = matrix[4][1] =-1
        matrix[3][4] = matrix[4][3] = 1


        F = [0] * 5

        F[0] = -(D[0][0])
        F[1] = -(D[0][1] - L[0])
        F[2] = -(D[1][0])
        F[3] = -(D[1][1] + L[0])

        F[4] = -(self.Points[1].aCoord[1]+D[1][1] - self.Points[0].aCoord[1] - D[0][1])

        return matrix, F
