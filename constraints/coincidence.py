from constraints.constraint import Constraint


class Coincidence(Constraint):
    """Совпадение 2-х точек"""
    def __init__(self, points):
        self.Points = points
        self.Ls = 2
        self.Name = "Совпадение"

    def get_description(self):
        return "Совпадение точек " + str(self.Points[0].v_return()) + " и " + str(self.Points[1].v_return())

    def LocalCon(self, deltas, lambdas):
        k = 6
        matrix = [0] * k
        for i in range(k):
            matrix[i] = [0] * k

        # Вклад в матрицу
        for i in range(4):
            matrix[i][i] = 1

        matrix[0][4] = matrix[4][0] = -1
        matrix[1][5] = matrix[5][1] = -1
        matrix[2][4] = matrix[4][2] = 1
        matrix[3][5] = matrix[5][3] = 1

        # Вклад в вектор F
        F = [0] * k

        F[0] = -(deltas[0][0] - lambdas[0])
        F[1] = -(deltas[0][1] - lambdas[1])
        F[2] = -(deltas[1][0] + lambdas[0])
        F[3] = -(deltas[1][1] + lambdas[1])

        F[4] = -(self.Points[1].aCoord[0] + deltas[1][0] - self.Points[0].aCoord[0] - deltas[0][0])
        F[5] = -(self.Points[1].aCoord[1] + deltas[1][1] - self.Points[0].aCoord[1] - deltas[0][1])

        return matrix, F
