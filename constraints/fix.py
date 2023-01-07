from constraints.constraint import Constraint


class Fixed(Constraint):
    """Класс фиксации"""
    def __init__(self, point):
        self.Points = point
        self.Ls = 2
        self.Name = "Фиксация"
        # Fixedlist.append(Points[0])
        # self.Flist = Fixedlist

    def get_description(self):
        return "Фиксация точки " + str(self.Points[0].v_return())

    def LocalCon(self, D, L):
        k = 4
        matrix = [0] * k
        for i in range(k):
            matrix[i] = [0] * k

        for i in range(2):
            matrix[i][i] = 1

        matrix[2][0] = matrix[0][2] = 0
        matrix[3][1] = matrix[1][3] = 0

        F = [0] * k

        F[0] = -(D[0][0] + L[0])
        F[1] = -(D[0][1] + L[1])
        F[2] = -(D[0][0])
        F[3] = -(D[0][1])

        return matrix, F

    def delConstrain(self):
        #self.Flist.remove(self.Points[0])
        #class_name = self.__class__.__name__
        #del self.Coords
        pass
