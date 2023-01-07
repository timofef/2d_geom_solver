from constraints.constraint import Constraint


class Fixed(Constraint):
    """Класс фиксации"""
    def __init__(self, point):
        self.Points = point
        self.x = point[0].aCoord[0]
        self.y = point[0].aCoord[1]
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

        matrix[0][2] = matrix[2][0] = 1
        matrix[1][3] = matrix[3][1] = 1

        F = [0] * k

        F[0] = -(D[0][0] + L[0])
        F[1] = -(D[0][1] + L[1])
        F[2] = -(self.Points[0].aCoord[0] + D[0][0] - self.x)
        F[3] = -(self.Points[0].aCoord[1] + D[0][1] - self.y)

        return matrix, F

    def delConstrain(self):
        #self.Flist.remove(self.Points[0])
        #class_name = self.__class__.__name__
        #del self.Coords
        pass
