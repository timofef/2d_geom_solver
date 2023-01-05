from constraints.constraint import Constraint


class Vertical(Constraint):
    """Класс вертикальности"""
    def __init__(self, Points):
        """Constructor"""
        self.Points = Points
        self.Ls = 1
        self.Name = "Вертикальность"

    def get_description(self):
        return "Вертикальность отрезка {" + str(self.Points[0].v_return()) + "; " + str(self.Points[1].v_return()) + "}"

    def LocalCon(self, deltas, lambdas):
        a = 5
        matrix = [0] * a
        for i in range(a):
            matrix[i] = [0] * a

        for i in range(4):
            matrix[i][i] = 1

        matrix[0][4] = matrix[4][0] = -1
        matrix[2][4] = matrix[4][2] = 1

        F = [0] * 5

        F[0] = -(deltas[0][0] - lambdas[0])
        F[1] = -(deltas[0][1])
        F[2] = -(deltas[1][0] + lambdas[0])
        F[3] = -(deltas[1][1])
        F[4] = -(self.Points[1].aCoord[0] + deltas[1][0] - self.Points[0].aCoord[0] - deltas[0][0])

        return matrix, F
