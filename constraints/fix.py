from constraints.constraint import Constraint


class Fixed(Constraint):
    """Класс фиксации"""
    def __init__(self, Points, Fixedlist):
        """Constructor"""
        self.Points = Points
        self.Ls = 0
        self.Name = "Фиксация"
        Fixedlist.append(Points[0])
        self.Flist = Fixedlist

    '''def LocalCon(self,D,L):
        a=2
        matrix=[0]*a
        for i in range(a):
            matrix[i]=[0]*a
        
        F=[0]*a

        F[0]=0
        F[1]=0

        return(matrix, F)'''

    def delConstrain(self):
        self.Flist.remove(self.Points[0])
        #class_name = self.__class__.__name__
        #del self.Coords
        pass
