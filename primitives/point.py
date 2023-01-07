from primitives.classes import Figure


class Point(Figure):
    """Класс точки"""

    def __init__(self, aCoord):
        self.aCoord = aCoord

    def v_return(self):
        return [round(self.aCoord[0], 10), round(self.aCoord[1], 10)]

    def delPoint(self):
        del self.aCoord
