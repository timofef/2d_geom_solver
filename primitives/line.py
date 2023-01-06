from primitives.classes import Figure


class Line(Figure):
    """Класс отрезка"""
    def __init__(self, point1, point2):
        self.oPoint1 = point1
        self.oPoint2 = point2

    def v_return(self):
        print('Point 1: \n\tX: {}; Y: {};\nPoint 2: \n\tX: {}; Y: {}.'.format(self.oPoint1.aCoord[0], self.oPoint1.aCoord[1], self.oPoint2.aCoord[0], self.oPoint2.aCoord[1]))

    def xy_return(self):
        aXcoord = [round(self.oPoint1.aCoord[0], 10), round(self.oPoint2.aCoord[0], 10)]
        aYcoord = [round(self.oPoint1.aCoord[1], 10), round(self.oPoint2.aCoord[1], 10)]
        return aXcoord, aYcoord

    def delLine(self):
        class_name = self.__class__.__name__
        #print('{} уничтожен'.format(class_name))
        del self.oPoint1
        del self.oPoint2
