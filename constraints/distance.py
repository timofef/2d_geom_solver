from constraints.constraint import Constraint


class Distance(Constraint):
    """Класс расстояния"""
    def __init__(self, points, d):
        """Constructor"""
        self.Points = points
        self.Ls = 1
        self.d = d
        self.Name = "Расстояние"

    def get_description(self):
        return "Расстояние между точками " + str(self.Points[0].v_return()) + " и " + str(self.Points[1].v_return())

    def LocalCon(self, deltas, lambdas):
        k = 5
        matrix = [0] * k
        for i in range(k):
            matrix[i] = [0] * k

        # Локальная матрица
        a = self.Points[1].aCoord[0] + deltas[1][0] - self.Points[0].aCoord[0] - deltas[0][0]
        b = self.Points[1].aCoord[1] + deltas[1][1] - self.Points[0].aCoord[1] - deltas[0][1]
        lam = lambdas[0]

        for i in range(4):
            matrix[i][i] = 1

        matrix[0][4] = matrix[4][0] = -2 * a
        matrix[1][4] = matrix[4][1] = -2 * b
        matrix[2][4] = matrix[4][2] = 2 * a
        matrix[3][4] = matrix[4][3] = 2 * b

        # Вклад в вектор F
        F = [0] * k

        F[0] = -(deltas[0][0] - 2 * lam * a)
        F[1] = -(deltas[0][1] - 2 * lam * b)
        F[2] = -(deltas[1][0] + 2 * lam * a)
        F[3] = -(deltas[1][1] + 2 * lam * b)

        F[4] = -(a*a + b*b - self.d*self.d)

        return matrix, F
