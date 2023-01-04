import sys
import pylab
from matplotlib.widgets import Button
from matplotlib.widgets import RadioButtons
from matplotlib.widgets import TextBox

import constraints.angle
import constraints.coincidence
import constraints.distance
import constraints.fix
import constraints.horizontal
import constraints.parallel
import constraints.perpendicular
import constraints.pointOnLine
import constraints.vertical
import objects.line
import objects.point
from objects import classes
from constraints import constraint
from PIL import Image
import numpy as np


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
FLAG_POL = 0  #
FLAG_CON = 0  # Задание совпадения точек
FLAG_ANG  = 0  # Задание угла
FLAG_PAR = 0  # Задание параллельности
FLAG_PER = 0  # Задание перпендикулярности


# Флаги режимов создания/удаления примитивов
FLAG_DEL = 0 # флаг для удаления точки
FLAG_CRE = 0 # флаг создания 0 - точка, 1 - прямая, 2 - прямая между существующими точками

PointCount = 0
PointInd = 0

cPointlist = []

global_point_list = []  # список точек
global_line_list = []  # список отрезков

XAXES = []  # размер рабочего поля по Х
YAXES = []  # размер рабочего поля по У

FLAG_DELCON = 0

Constraints = []
Fixedlist = []

EPS = 0.01


def submit(expression):
    """
    Update the plotted function to the new math *expression*.

    *expression* is a string using "t" as its independent variable, e.g.
    "t ** 3".
    """
    ydata = eval(expression)
    return ydata


def check_time(curtime, mode):
    global CLICK_LAST_TIME, PICK_LAST_TIME
    if mode == 0: # Режим просто для клика
        delta = curtime - CLICK_LAST_TIME
        if delta.microseconds < 100000:
            #print(delta)
            return 0
        else:
            CLICK_LAST_TIME = curtime
            return 1
    if mode == 1:  # Режим для PICK
        delta = curtime - PICK_LAST_TIME
        if delta.microseconds < 100000:
            #print(delta)
            return 0
        else:
            PICK_LAST_TIME = curtime
            return 1


# Сброс флагов режима
def reset_flags():
    global FLAG_HOR, FLAG_VER, FLAG_FIX, FLAG_DIS, FLAG_POL, FLAG_CON, FLAG_ANG, FLAG_PAR, FLAG_PER
    FLAG_HOR = 0
    FLAG_VER = 0
    FLAG_FIX = 0
    FLAG_DIS = 0
    FLAG_POL = 0
    FLAG_CON = 0
    FLAG_ANG  = 0
    FLAG_PAR = 0
    FLAG_PER = 0


# Итерации метода Ньютона
def Redoer():
    global ULTRACOUNTER
    if len(Constraints) != 0:
        tmp = 2 * len(global_point_list)
        for con in Constraints:
            tmp += con.getLs()
        deltas = [0] * tmp  #starting deltas

        ULTRACOUNTER = 0

        while True:
            matrix, F = assemble_slae(deltas)
            Adders, flag = solve_slae(matrix, F)
            if flag:
                print("Не могу пересчитать, отменяюсь")
                return 1
            if abs(max(Adders, key=abs)) < EPS:
                for i in range(len(global_point_list)):
                    global_point_list[i].aCoord[0] += deltas[i * 2]
                    global_point_list[i].aCoord[1] += deltas[i * 2 + 1]
                break

            for i in range(len(deltas)):
                deltas[i] += Adders[i]
    return 0


#  Формирование матрицы системы
def assemble_slae(deltas):
    d_len = len(deltas)
    matrix = [0] * d_len
    for i in range(d_len):
        matrix[i] = [0] * d_len
    f = [0] * d_len

    # Заносим координаты всех точек в список
    coordinates = []
    total_coordinates = 2 * len(global_point_list)
    for i in range(0, total_coordinates, 2):
        coordinates.append([deltas[i], deltas[i + 1]])

    # Формируем список лямбд для всех имеющихся ограничений
    lamdas = []
    for con in Constraints:
        tmp = [0] * con.getLs()
        for i in range(con.getLs()):
            tmp[i] = deltas[total_coordinates]
            total_coordinates += 1
        lamdas.append(tmp)

    lamda_shift = 0
    for con in Constraints:
        if con.Ls != 0:
            D = []
            L = lamdas[Constraints.index(con)]
            T = []  # Координаты в глобальной матрице

            for cp in con.Points:
                for gp in global_point_list:
                    if cp == gp:
                        D.append(coordinates[global_point_list.index(gp)])
                        T.append(2 * global_point_list.index(gp))
                        T.append(2 * global_point_list.index(gp) + 1)
                        break
            for i in range(len(L)):
                T.append(2 * len(global_point_list) + lamda_shift + i)
            l_matrix, l_F = con.LocalCon(D, L)
            lamda_shift += len(L)

            # Ансамблирование
            for i in range(len(T)):
                for j in range(len(T)):
                    matrix[T[i]][T[j]] += l_matrix[i][j]
                f[T[i]] += l_F[i]
    
    for i in range(2 * len(global_point_list)):  # для точек, не попавших в ограничения
        if matrix[i][i] == 0 and f[i] == 0:
            matrix[i][i] = 1

    for i in range(len(global_point_list)):
        if global_point_list[i] in Fixedlist:
            for j in range(len(deltas)):
                matrix[i*2][j] = 0
                matrix[i*2 + 1][j] = 0
                matrix[j][i*2] = 0
                matrix[j][i*2 + 1] = 0
            matrix[i*2][i*2] = 1
            matrix[i*2 + 1][i*2 + 1] = 1
            f[i*2] = 0
            f[i*2 + 1] = 0

    return matrix, f


# Решение СЛАУ методом Гаусса
def solve_slae(matrix, f):
    # Решение СЛАУ методом Гаусса
    global ULTRACOUNTER
    alen = len(f)
    flag = 0
    for h in range(alen - 1):
        for j in range(h + 1, alen):
            if matrix[j][h] != 0:
                if -0.001 < matrix[h][h] < 0.001:
                    flag = 1
                    break
                else:
                    m = matrix[j][h] / matrix[h][h]
                try:
                    for i in range(h, alen):
                        matrix[j][i] -= m * matrix[h][i]
                    try:
                        f[j] -= m * f[h]
                        flag = 0
                    except:
                        flag = 1
                    flag = 0
                except:
                    flag = 1
    result = [0] * alen

    if flag:
        return result, flag

    if -0.001 < matrix[alen - 1][alen - 1] < 0.001:
        return result, flag

    for h in range(alen - 1, -1, -1):
        m = 0
        if not(-0.001 < matrix[h][h] < 0.001):
            for i in range(h + 1, alen):
                m += matrix[h][i] * result[i]
            result[h] = (f[h] - m) / matrix[h][h]
        else:
            # print("Что-то определенно не так, скорее всего - переопределенность, а именно - излишняя фиксация")
            ULTRACOUNTER += 1
            if ULTRACOUNTER > 9:
                flag = 1

    return result, flag


# Отрисовка эскиза
def draw_plot(): # Перерисовка поля graph_axes
    global global_point_list, global_line_list, cPointlist

    cPointlist = []
    cPointlist = global_point_list.copy()
    if Redoer():
        global_point_list = cPointlist.copy()

    aX = []
    aY = []
    graph_axes.clear()
    graph_axes.set_xlim(XAXES)
    graph_axes.set_ylim(YAXES)

    for oPoint in global_point_list:
        aX.append(oPoint.v_return()[0])
        aY.append(oPoint.v_return()[1])

    for oLine in global_line_list:
        x, y = oLine.xy_return()
        graph_axes.plot(x, y, marker = 'o', picker = True, pickradius = 5)

    graph_axes.plot(aX, aY, 'o', color='red', picker = True, pickradius = 5, zorder = 2.5)
    graph_axes.grid()
    pylab.draw()


# Обработка кликов в рабочем поле
def on_click(event):
    """Нажатие по пустой области"""
    global FLAG_CRE, PointCount, global_point_list
    if check_time(datetime.now(), 0):
        if event.button == 3 and FLAG_CRE == 0:  # Создание точки
            add_point(event)
        elif event.button == 3 and FLAG_CRE == 1:  # Создание отрезка
            if PointCount == 0:  # Выбрана первая точка
                add_point(event)
                PointCount = 1
            elif PointCount == 1: # Выбрана вторая точка
                add_point(event)
                add_line(global_point_list[-2], global_point_list[-1])
                PointCount = 0


def on_pick(event):
    """Обработчик клика по существующему примитиву"""
    global tmp, global_point_list, graph_axes, FLAG_DEL, FLAG_CRE, PointCount, PointInd, FLAG_FIX, \
        Constraints, FLAG_HOR, FLAG_VER, FLAG_DIS, FLAG_POL, FLAG_CON, FLAG_ANG, FLAG_PAR, FLAG_PER
    if check_time(datetime.now(), 1):
        thisline = event.artist
        xdata = thisline.get_xdata()
        ydata = thisline.get_ydata()
        ind = event.ind

        msx, msy = event.mouseevent.xdata, event.mouseevent.ydata
        dist = np.sqrt((np.array(xdata)-msx)**2+(np.array(ydata)-msy)**2)
        mindist = min(dist)

        dist_kostyl = 1.0  # Поиск индекса ближайшей точки
        pointInd = -1
        for oPoint in global_point_list:
            dist_kostyl = np.sqrt((oPoint.aCoord[0]-msx)**2+(oPoint.aCoord[1]-msy)**2)
            if dist_kostyl < 0.1:
                pointInd = global_point_list.index(oPoint)
        dist_kostyl = 1.0

        # Если кликнули на отрезок
        if mindist > 0.1 and pointInd == -1:
            if FLAG_DEL:  # при удалении линии
                for oLine in global_line_list:
                    if get_line_by_coord(oLine, xdata, ydata):
                        delete_line(oLine)
                        FLAG_DEL = 0
                        button_del.color = 'white'

            if FLAG_HOR: # горизонтальность
                for oLine in global_line_list:
                    if get_line_by_coord(oLine, xdata, ydata):
                        tmp = [oLine.oPoint1, oLine.oPoint2]
                        tmp2= constraints.horizontal.Horizontal(tmp)
                        Constraints.append(tmp2)
                        draw_plot()
                        FLAG_HOR = 0

            if FLAG_VER: # вертикальность
                for oLine in global_line_list:
                    if get_line_by_coord(oLine, xdata, ydata):
                        tmp = [oLine.oPoint1, oLine.oPoint2]
                        tmp2= constraints.vertical.Vertical(tmp)
                        Constraints.append(tmp2)
                        draw_plot()
                        FLAG_VER = 0
            

            if FLAG_PAR == 2:  # параллельность
                for oLine in global_line_list:
                    if get_line_by_coord(oLine, xdata, ydata):
                        tmp.append(oLine.oPoint1)
                        tmp.append(oLine.oPoint2)
                        tmp2 = constraints.parallel.Parallel(tmp)
                        Constraints.append(tmp2)
                        draw_plot()
                        FLAG_PAR = 0
            if FLAG_PAR == 1:
                for oLine in global_line_list:
                    if get_line_by_coord(oLine, xdata, ydata):
                        tmp = [oLine.oPoint1, oLine.oPoint2]
                        print("Первая линия выбрана, ждем вторую")
                        FLAG_PAR = 2

            
            if FLAG_PER == 2:  # перпендикулярность
                for oLine in global_line_list:
                    if get_line_by_coord(oLine, xdata, ydata):
                        tmp.append(oLine.oPoint1)
                        tmp.append(oLine.oPoint2)
                        tmp2 = constraints.perpendicular.Perpendicular(tmp)
                        Constraints.append(tmp2)
                        draw_plot()
                        FLAG_PER = 0
            if FLAG_PER == 1:
                for oLine in global_line_list:
                    if get_line_by_coord(oLine, xdata, ydata):
                        tmp = [oLine.oPoint1, oLine.oPoint2]
                        print("Первая линия выбрана, ждем вторую")
                        FLAG_PER = 2

            if FLAG_POL == 1: # точка на линии
                for oLine in global_line_list:
                    if get_line_by_coord(oLine, xdata, ydata):
                        tmp = [oLine.oPoint1, oLine.oPoint2]
                        print("Линия выбрана, теперь точка")
                        FLAG_POL = 2

            if FLAG_ANG == 2:  # угол, выбор второго отрезка и задание угла
                for oLine in global_line_list:
                    if get_line_by_coord(oLine, xdata, ydata):
                        tmp.append(oLine.oPoint1)
                        tmp.append(oLine.oPoint2)
                        print("Отлично, теперь угол")
                        u = input()
                        flag = 1
                        try:
                            u = float(u)
                        except:
                            flag = 0
                        if flag:
                            tmp2 = constraints.angle.Angle(tmp, u)
                            Constraints.append(tmp2)
                            draw_plot()
                            FLAG_ANG = 0
                        else:
                            print("Плохо введен угол, отменяемся")
                            FLAG_ANG = 0
            if FLAG_ANG == 1:
                for oLine in global_line_list:
                    if get_line_by_coord(oLine, xdata, ydata):
                        tmp = [oLine.oPoint1, oLine.oPoint2]
                        print("Первая линия выбрана, ждем вторую")
                        FLAG_ANG = 2
        # Если кликнули на точку
        else:
            points = tuple(zip(xdata[ind], ydata[ind]))
            if FLAG_DEL: # при удалении точки
                if pointInd != -1:
                    delete_point(global_point_list[pointInd])
                    FLAG_DEL = 0
                    button_del.color = 'white'
                    pointInd = -1
                else:
                    for oPoint in global_point_list:
                        if oPoint.v_return() == [round(points[0][0],10),round(points[0][1],10)]:
                            delete_point(oPoint)
                            FLAG_DEL = 0
                            button_del.color = 'white'

            if FLAG_CRE == 2 and event.mouseevent.button == 3:  # Добавление отрезка по 2-м существующим точкам
                if PointCount == 0:
                    if pointInd != -1:
                        PointInd = global_point_list.index(global_point_list[pointInd])
                        print("Нашли первую точку")
                        PointCount = 1
                        pointInd = -1
                        return 0
                    else:
                        for oPoint in global_point_list:
                            if oPoint.v_return() == [round(points[0][0], 10), round(points[0][1], 10)]:
                                PointInd = global_point_list.index(oPoint)
                                print("Нашли первую точку")
                                PointCount = 1
                                return 0
                if PointCount == 1:
                    if pointInd != -1:
                        if global_point_list[PointInd] != global_point_list[pointInd]:
                            print("Нашли вторую точку")
                            add_line(global_point_list[PointInd], global_point_list[pointInd])
                            PointCount = 0
                            PointInd = 0
                            pointInd = -1
                            return 0
                    else:
                        for oPoint in global_point_list:
                            if oPoint.v_return() == [round(points[0][0], 10), round(points[0][1], 10)]:
                                if global_point_list[PointInd] != oPoint:
                                    print ("Нашли вторую точку")
                                    add_line (global_point_list[PointInd], oPoint)
                                    PointCount = 0
                                    PointInd = 0
                                    return 0
            if FLAG_FIX == 1:
                if pointInd != -1:
                    tmp = [global_point_list[pointInd]]
                    tmp2 = constraints.fix.Fixed(tmp, Fixedlist)
                    Constraints.append(tmp2)
                    draw_plot()
                    print("Точка зафиксирована")
                    FLAG_FIX = 0
                    pointInd = -1
                else:
                    for oPoint in global_point_list:
                        if oPoint.v_return() == [round(points[0][0], 10), round(points[0][1], 10)]:
                            tmp = [oPoint]
                            tmp2 = constraints.fix.Fixed(tmp, Fixedlist)
                            Constraints.append(tmp2)
                            draw_plot()
                            print("Точка зафиксирована")
                            FLAG_FIX = 0
            
            if FLAG_DIS == 2:
                if pointInd != -1:
                    tmp.append(global_point_list[pointInd])
                    print("Вторая точка выбрана, введите расстояние")
                    # axbox1 = fig.add_axes([0.07, 0.63, 0.9, 0.08])
                    # text_box1 = TextBox(axbox1, '', textalignment="center")
                    # text_box1.set_val("Введите данные")  # Trigger `submit` with the initial string.
                    # d = text_box1.on_submit(submit)
                    # print(d)
                    d = input()
                    flag = 1
                    try:
                        d = float(d)
                    except:
                        flag = 0
                    if flag:
                        tmp2= constraints.distance.Distance(tmp, d)
                        Constraints.append(tmp2)
                        draw_plot()
                        pointInd = -1
                        FLAG_DIS = 0
                    else:
                        print("Плохо введено расстояние, отменяемся")
                        pointInd = -1
                        FLAG_DIS = 0
                else:
                    for oPoint in global_point_list:
                        if oPoint.v_return() == [round(points[0][0], 10), round(points[0][1], 10)]:
                            tmp.append(oPoint)
                            print("Вторая точка выбрана, введите расстояние")
                            d = input()
                            flag=1

                            try:
                                d = float(d)
                            except:
                                flag = 0

                            if flag:
                                tmp2 = constraints.distance.Distance(tmp, d)
                                Constraints.append(tmp2)
                                draw_plot()
                                FLAG_DIS = 0
                            else:
                                print("Плохо введено расстояние, отменяемся")
                                FLAG_DIS = 0

            if FLAG_DIS == 1:
                if pointInd != -1:
                    tmp = []
                    tmp.append(global_point_list[pointInd])
                    print("Первая точка выбрана, ждем вторую")
                    FLAG_DIS = 2
                    pointInd = -1
                else:
                    for oPoint in global_point_list:
                        if oPoint.v_return() == [round(points[0][0], 10), round(points[0][1], 10)]:
                            tmp = []
                            tmp.append(oPoint)
                            print("Первая точка выбрана, ждем вторую")
                            FLAG_DIS = 2
            

            if FLAG_CON == 2:
                if pointInd != -1:
                    tmp.append(global_point_list[pointInd])
                    #print("Вторая точка выбрана, введите расстояние")
                    tmp2= constraints.coincidence.Coincidence(tmp)
                    Constraints.append(tmp2)
                    draw_plot()
                    FLAG_CON = 0
                    pointInd = -1
                else:
                    for oPoint in global_point_list:
                        if oPoint.v_return() == [round(points[0][0], 10), round(points[0][1], 10)]:
                            tmp.append(oPoint)
                            #print("Вторая точка выбрана, введите расстояние")
                            tmp2= constraints.coincidence.Coincidence(tmp)
                            Constraints.append(tmp2)
                            draw_plot()
                            FLAG_CON = 0


            if FLAG_CON == 1:
                if pointInd != -1:
                    tmp = []
                    tmp.append(global_point_list[pointInd])
                    print("Первая точка выбрана, ждем вторую")
                    FLAG_CON = 2
                    pointInd = -1
                else:
                    for oPoint in global_point_list:
                        if oPoint.v_return() == [round(points[0][0],10),round(points[0][1],10)]:
                            tmp = []
                            tmp.append(oPoint)
                            print("Первая точка выбрана, ждем вторую")
                            FLAG_CON = 2


            if FLAG_POL == 2:
                if pointInd != -1:
                    tmp.append(global_point_list[pointInd])
                    #print("Вторая точка выбрана, введите расстояние")
                    tmp2= constraints.pointOnLine.PointOnLine(tmp)
                    Constraints.append(tmp2)
                    draw_plot()
                    FLAG_POL = 0
                    pointInd = -1
                else:
                    for oPoint in global_point_list:
                        if oPoint.v_return() == [round(points[0][0], 10), round(points[0][1], 10)]:
                            tmp.append(oPoint)
                            #print("Вторая точка выбрана, введите расстояние")
                            tmp2 = constraints.pointOnLine.PointOnLine(tmp)
                            Constraints.append(tmp2)
                            draw_plot()
                            FLAG_POL = 0


def get_line_by_coord(line, x, y):
    a = line.xy_return()[0][0] == round(x[0], 10)
    b = line.xy_return()[0][1] == round(x[1], 10)
    c = line.xy_return()[1][0] == round(y[0], 10)
    d = line.xy_return()[1][1] == round(y[1], 10)
    return a and b and c and d

# Создание примитивов
def add_point(event):
    '''Добавление новой точки'''
    global global_point_list
    new_point = objects.point.Point([event.xdata, event.ydata])
    global_point_list.append(new_point)
    draw_plot()


def add_line(p1, p2):
    '''Добавление нового отрезка'''
    global global_line_list
    new_line = objects.line.Line(p1, p2)
    global_line_list.append(new_line)
    draw_plot()


# Удаление примитивов
def delete_point(oPoint):
    '''Удаление точки'''
    global global_point_list, global_line_list
    # Если к точке привязаны прямые, то они тоже стираются
    aLinelist_copy = []
    for oLine in global_line_list:
        print(oLine.xy_return())
        if oPoint == oLine.oPoint1 or oPoint == oLine.oPoint2:
            aLinelist_copy.append(oLine)
    for oLine1 in aLinelist_copy:
        global_line_list.remove(oLine1)
        oLine1.delLine()

    # конец блока удаления связанных прямых

    while True:
        ttt = 0
        for con in Constraints:
            if oPoint in con.Points:
                ttt = 1
                Constraints.remove(con)
        if ttt == 0:
            break
    print(len(Constraints))

    global_point_list.remove(oPoint)
    oPoint.delPoint()
    draw_plot()


def delete_line(oLine):
    global global_point_list, global_line_list
    
    oPoint1 = global_point_list[global_point_list.index(oLine.oPoint1)]
    oPoint2 = global_point_list[global_point_list.index(oLine.oPoint2)]
    iP1cnt = -1
    iP2cnt = -1
    for oLine_ in global_line_list: # Считаем кол-во точек из Линии в списке линий
        if oPoint1 == oLine_.oPoint1 or oPoint1 == oLine_.oPoint2:
            iP1cnt +=1
        if oPoint2 == oLine_.oPoint1 or oPoint2 == oLine_.oPoint2:
            iP2cnt += 1
    global_line_list.remove(oLine)
    oLine.delLine()
    # Если точка Линии больше нигде не используется, то её стираем, иначе не трогаем
    if not iP1cnt:
        while True:
            ttt = 0
            for con in Constraints:
                if oPoint1 in con.Points:
                    ttt = 1
                    Constraints.remove(con)
            if ttt == 0:
                break
        print(len(Constraints))

        global_point_list.remove(oPoint1)
        oPoint1.delPoint()
    if not iP2cnt:
        while True:
            ttt = 0
            for con in Constraints:
                if oPoint2 in con.Points:
                    ttt = 1
                    Constraints.remove(con)
            if ttt == 0:
                break
        print(len(Constraints))

        global_point_list.remove(oPoint2)
        oPoint2.delPoint()
    draw_plot()


# Переключение режимов создания/удаления примитивов
def choose_delete_primitive_mode(event):
    """Обработчик события для кнопки "Удалить точку\""""
    global FLAG_DEL
    if FLAG_DEL == 1:
        FLAG_DEL = 0
        button_del.color = 'white'
    else:
        FLAG_DEL = 1
        button_del.color = 'grey'


def choose_delete_constraint_mode(event):
    """Обработчик события для кнопки "Удалить ограничение\""""
    global FLAG_DELCON
    if FLAG_DELCON == 1:
        FLAG_DELCON = 0
        button_delcon.color = 'white'
    else:
        if len(Constraints) != 0:
            FLAG_DELCON = 1
            button_delcon.color = 'grey'
            print("Список ограничений:")
            i = 0
            for con in Constraints:
                tmp = []
                for point in con.Points:
                    tmp.append(point.v_return())
                print(i, ":", con.Name, "между точками:", tmp)
                i += 1
            print("Ожидаю, какое по номеру ограничение удалить")
            flag = 1
            todel = input()
            try:
                todel = int(todel)
            except:
                flag = 0
            if flag:
                try:
                    Constraints[todel].delConstrain()
                    Constraints.remove(Constraints[todel])
                except:
                    print("Ограничения с таким индексом нет")
                print("Новый список ограничений:")
                i = 0
                for con in Constraints:
                    tmp = []
                    for point in con.Points:
                        tmp.append(point.v_return())
                    print(i,":",con.Name,"между точками:",tmp)
                    i += 1
            FLAG_DELCON = 0
            button_delcon.color = 'white'
        else:
            print("Ограничений нет, удалять нечего")


def choose_creation_mode(label):
    """Выбор примитива для ввода"""
    global FLAG_CRE
    if radiobuttons_creationtype.value_selected == 'Точка':
        FLAG_CRE = 0
    elif radiobuttons_creationtype.value_selected == 'Прямая':
        FLAG_CRE = 1
    else:
        FLAG_CRE = 2


def onButtonFpClicked(event):
    """Обработка нажатия на кнопки ограничений"""
    global FLAG_FIX

    if FLAG_FIX == 1:  # Начать заново при повторном клике на кнопку
        FLAG_FIX = 0
    else:
        text_box.set_val("Выберите точку для фиксации")
        reset_flags()
        FLAG_FIX = 1


def onButtonDisClicked(event):
    """Наложение расстояния"""
    global FLAG_DIS

    if FLAG_DIS == 1 or FLAG_DIS == 2:
        FLAG_DIS = 0
    else:
        text_box.set_val("Выберите точки для наложения расстояния")
        reset_flags()
        FLAG_DIS = 1


def onButtonPolClicked(event):
    """Наложение точки на линии"""
    global FLAG_POL

    if FLAG_POL == 1 or FLAG_POL == 2:
        FLAG_POL = 0
    else:
        text_box.set_val("Выберите линию и точку для принадлежности")
        reset_flags()
        FLAG_POL = 1


def onButtonConClicked(event):
    """Наложение совпадения точек"""
    global FLAG_CON

    if FLAG_CON == 1 or FLAG_CON == 2:
        FLAG_CON = 0
    else:
        text_box.set_val("Выберите точки для совпадения")
        reset_flags()
        FLAG_CON = 1


def onButtonAnClicked(event):
    """Задание угла"""
    global FLAG_ANG

    if FLAG_ANG == 1 or FLAG_ANG == 2:
        FLAG_ANG = 0
    else:
        text_box.set_val("Выберите линии для наложения угла")
        reset_flags()
        FLAG_ANG = 1


def onButtonParClicked(event):
    """Наложение параллельности"""
    global FLAG_PAR

    if FLAG_PAR == 1 or FLAG_PAR == 2:
        FLAG_PAR = 0
    else:
        text_box.set_val("Выберите прямые для параллельности")
        reset_flags()
        FLAG_PAR = 1


def onButtonPerClicked(event):
    """Наложение перпендикулярности"""
    global FLAG_PER

    if FLAG_PER == 1 or FLAG_PER == 2:
        FLAG_PER = 0
    else:
        text_box.set_val("Выберите прямые для перпендикулярности")
        reset_flags()
        FLAG_PER = 1


def onButtonVerClicked(event):
    global FLAG_VER
    
    if FLAG_VER == 1:
        FLAG_VER = 0
    else:
        text_box.set_val("Выберите линию наложения вертикальности")
        reset_flags()
        FLAG_VER = 1


def onButtonHorClicked(event):
    global FLAG_HOR

    if FLAG_HOR == 1:
        FLAG_HOR = 0
    else:
        text_box.set_val("Выберите линию наложения горизонтальности")
        reset_flags()
        FLAG_HOR = 1


if __name__ == "__main__":
    XAXES = [int(sys.argv[1]), int(sys.argv[2])]  # Command line arguments for field size
    YAXES = [int(sys.argv[3]), int(sys.argv[4])]
    fig, graph_axes = pylab.subplots()

    graph_axes.set_xlim(XAXES)
    graph_axes.set_ylim(YAXES)
    graph_axes.aspect = 1
    graph_axes.grid()
    pylab.tight_layout()

    # Размеры окна и место для кнопок
    fig.set_size_inches(16, 9)
    fig.subplots_adjust(left=0.07, right=0.8, top=0.95, bottom=0.15)

    # Кнопки для ограничений
    axes_button_fp = pylab.axes([0.81, 0.865, 0.085, 0.085])
    img_fp = Image.open("./icons/fixpoint.png")
    button_fp = Button(axes_button_fp, None, img_fp)
    button_fp.on_clicked(onButtonFpClicked)

    axes_button_ver = pylab.axes([0.905, 0.865, 0.085, 0.085])
    img_ver = Image.open("./icons/verticality.png")
    button_ver = Button(axes_button_ver, None, img_ver)
    button_ver.on_clicked(onButtonVerClicked)

    axes_button_dis = pylab.axes([0.81, 0.77, 0.085, 0.085])
    img_dis = Image.open("./icons/distance.png")
    button_dis = Button(axes_button_dis, None, img_dis)
    button_dis.on_clicked(onButtonDisClicked)

    axes_button_hor = pylab.axes([0.905, 0.77, 0.085, 0.085])
    img_hor = Image.open("./icons/horizontality.png")
    button_hor = Button(axes_button_hor, None, img_hor)
    button_hor.on_clicked(onButtonHorClicked)

    axes_button_pol = pylab.axes([0.905, 0.675, 0.085, 0.085])
    img_pol = Image.open("./icons/pointonline.png")
    button_pol= Button(axes_button_pol, None, img_pol)
    button_pol.on_clicked(onButtonPolClicked)

    axes_button_con = pylab.axes([0.81, 0.675, 0.085, 0.085])
    img_con = Image.open("./icons/concidence.png")
    button_con= Button(axes_button_con, None, img_con)
    button_con.on_clicked(onButtonConClicked)

    axes_button_an = pylab.axes([0.81, 0.58, 0.085, 0.085])
    img_an = Image.open("./icons/angle.png")
    button_an= Button(axes_button_an, None, img_an)
    button_an.on_clicked(onButtonAnClicked)

    axes_button_par = pylab.axes([0.905, 0.58, 0.085, 0.085])
    img_par = Image.open("./icons/parallelism.png")
    button_par = Button(axes_button_par, None, img_par)
    button_par.on_clicked(onButtonParClicked)

    axes_button_per = pylab.axes([0.81, 0.485, 0.085, 0.085])
    img_per = Image.open("./icons/perpendicularity.png")
    button_per = Button(axes_button_per, None, img_per)
    button_per.on_clicked(onButtonPerClicked)

    # Кнопки для добаления примитивов
    axes_radiobuttons = pylab.axes([0.829, 0.375, 0.14, 0.1])
    radiobuttons_creationtype = RadioButtons(axes_radiobuttons, ['Точка', 'Прямая', 'Прямая между точек'])
    radiobuttons_creationtype.on_clicked(choose_creation_mode)

    # 0.39 0.295
    axes_button_del = pylab.axes([0.829, 0.285, 0.14, 0.08])
    button_del = Button(axes_button_del, 'Удалить примитив', color="white")
    button_del.on_clicked(choose_delete_primitive_mode)

    axes_button_delcon = pylab.axes([0.829, 0.195, 0.14, 0.08])
    button_delcon = Button(axes_button_delcon, 'Удалить ограничение', color="white")
    button_delcon.on_clicked(choose_delete_constraint_mode)

    # Обработчики кликов
    fig.canvas.mpl_connect('button_press_event', on_click)
    fig.canvas.mpl_connect('pick_event', on_pick)

    axbox = fig.add_axes([0.07, 0.03, 0.6, 0.08])
    text_box = TextBox(axbox, None, textalignment="center")
    text_box.set_val("Введите данные")  # Trigger `submit` with the initial string.

    pylab.show()