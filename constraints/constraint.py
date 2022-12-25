class Constraint:
    """Базовый класс условия"""
    def __init__(self, points):
        self.Points = points
        self.Ls = 0
        self.Name = "Base"

    def getLs(self):
        return self.Ls

    def delConstrain(self):  
        #class_name = self.__class__.__name__  
        #del self.Coords
        pass
