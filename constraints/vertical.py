from constraints.constraint import Constraint


class Vertical(Constraint):
    """Класс вертикальности"""
    def __init__(self, Points):
        """Constructor"""
        self.Points = Points
        self.Ls = 1
        self.Name = "Вертикальность"

    def LocalCon(self, deltas, lambdas):
        a = 5
        matrix = [0] * a
        for i in range(a):
            matrix[i] = [0] * a

        matrix[0][0] = 1
        matrix[1][1] = 1
        matrix[2][2] = 1
        matrix[3][3] = 1

        matrix[0][4] = -1
        matrix[2][4] = 1
        matrix[4][0] = -1
        matrix[4][2] = 1

        F=[0]*5

        F[0] = -(deltas[0][0] - lambdas[0])
        F[1] = -(deltas[0][1])
        F[2] = -(deltas[1][0] + lambdas[0])
        F[3] = -(deltas[1][1])
        F[4] = -(self.Points[1].aCoord[0] + deltas[1][0] - self.Points[0].aCoord[0] - deltas[0][0])

        return matrix, F
