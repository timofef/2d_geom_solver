from datetime import datetime

CLICK_LAST_TIME = datetime.now()
PICK_LAST_TIME = datetime.now()


tmp = []

ULTRACOUNTER = 0


# Флаги текущего режима
FLAG_HOR = 0  # Задание горизонтальности
FLAG_VER = 0  # Задание вертикальности
FLAG_FIX = 0  # Фиксация точки
FLAG_DIS = 0  # Расстояние
FLAG_POL = 0  # Точка на прямой
FLAG_CON = 0  # Задание совпадения точек
FLAG_ANG = 0  # Задание угла
FLAG_PAR = 0  # Задание параллельности
FLAG_PER = 0  # Задание перпендикулярности


# Флаги режимов создания/удаления примитивов
FLAG_DEL = 0  # флаг для удаления точки
FLAG_CRE = 0  # флаг создания 0 - точка, 1 - прямая, 2 - прямая между существующими точками

PointCount = 0
PointInd = 0

cPointlist = []

global_point_list = []  # список точек
global_line_list = []  # список отрезков

XAXES = [0, 10]  # размер рабочего поля по Х
YAXES = [0, 10]  # размер рабочего поля по У

FLAG_DELCON = 0

Constraints = []
Fixedlist = []

EPS = 0.001
