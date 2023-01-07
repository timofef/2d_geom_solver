import matplotlib.pyplot as plt
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
import primitives.line
import primitives.point
from PIL import Image
import numpy as np

from globals import *


# Задержка при нажатиях (для избежания случайных нажатий)
def check_time(curtime, mode):
    global CLICK_LAST_TIME, PICK_LAST_TIME
    if mode == 0:  # Режим просто для клика
        delta = curtime - CLICK_LAST_TIME
        if delta.microseconds < 100000:
            return 0
        else:
            CLICK_LAST_TIME = curtime
            return 1
    if mode == 1:  # Режим для PICK
        delta = curtime - PICK_LAST_TIME
        if delta.microseconds < 100000:
            return 0
        else:
            PICK_LAST_TIME = curtime
            return 1


# Сброс флагов режима
def reset_flags():
    global FLAG_HOR, FLAG_VER, FLAG_FIX, FLAG_DIS, \
        FLAG_POL, FLAG_CON, FLAG_ANG, FLAG_PAR, FLAG_PER
    FLAG_HOR = 0
    FLAG_VER = 0
    FLAG_FIX = 0
    FLAG_DIS = 0
    FLAG_POL = 0
    FLAG_CON = 0
    FLAG_ANG = 0
    FLAG_PAR = 0
    FLAG_PER = 0


# Перерасчёт координат точек
def update_primitives():
    global ULTRACOUNTER

    if len(Constraints) != 0:
        tmp = 2 * len(global_point_list)
        l = 0
        for con in Constraints:
            l += con.getLs()
        deltas = [0] * (tmp + l)

        if len(global_point_list) * 2 < l + len(Fixedlist) * 2:
            print("Переопределённость. Удалите последнее добавленное ограничение")
            return 1

        ULTRACOUNTER = 0

        # Итерации метода Ньютона
        while True:
            matrix, f = assemble_slae(deltas)
            new_deltas, flag = solve_slae(matrix, f)
            if flag:
                message_box.set_val('Не могу пересчитать, отменяюсь')
                # print("Не могу пересчитать, отменяюсь")
                return 1
            if abs(max(new_deltas, key=abs)) < EPS:
                for i in range(len(global_point_list)):
                    global_point_list[i].aCoord[0] += deltas[i * 2]
                    global_point_list[i].aCoord[1] += deltas[i * 2 + 1]
                break

            for i in range(len(deltas)):
                deltas[i] += new_deltas[i]
    return 0


#  Ансамблирование
def assemble_slae(deltas):
    d_len = len(deltas)
    matrix = [0] * d_len
    for i in range(d_len):
        matrix[i] = [0] * d_len
    f = [0] * d_len

    # Заносим координаты всех точек в список
    delta_coordinates = []
    total_coordinates = 2 * len(global_point_list)
    for i in range(0, total_coordinates, 2):
        delta_coordinates.append([deltas[i], deltas[i + 1]])

    # Формируем список лямбд для всех имеющихся ограничений
    delta_lamdas = []
    for con in Constraints:
        tmp = [0] * con.getLs()
        for i in range(con.getLs()):
            tmp[i] = deltas[total_coordinates]
            total_coordinates += 1
        delta_lamdas.append(tmp)

    lamda_shift = 0
    for con in Constraints:
        if con.Ls != 0:
            D = []  # Дельты для перерасчёта локальных матриц
            L = delta_lamdas[Constraints.index(con)]
            T = []  # Индекс в глобальной матрице

            # Получение дельт для точек, задействованных в данном ограничении
            for cp in con.Points:
                point_index = global_point_list.index(cp)
                D.append(delta_coordinates[point_index])
                T.append(2 * point_index)
                T.append(2 * point_index + 1)
            # Получение дельт лямбда для данного ограничения и сдвига в индексах
            # глобального вектора дельт для них
            for i in range(len(L)):
                T.append(2 * len(global_point_list) + lamda_shift + i)
            lamda_shift += len(L)

            # Получение локальной матрицы для данного ограничения
            local_matrix, local_f = con.LocalCon(D, L)

            # Ансамблирование
            for i in range(len(T)):
                for j in range(len(T)):
                    matrix[T[i]][T[j]] += local_matrix[i][j]
                f[T[i]] += local_f[i]

    # Точки без ограничений
    for i in range(2 * len(global_point_list)):
        if matrix[i][i] == 0 and f[i] == 0:
            matrix[i][i] = 1

    # Учёт фиксированных точек
    for i in range(len(global_point_list)):
        if global_point_list[i] in Fixedlist:
            for j in range(len(deltas)):
                matrix[i*2][j] = matrix[j][i*2] = 0
                matrix[i*2 + 1][j] = matrix[j][i*2 + 1] = 0
            matrix[i*2][i*2] = matrix[i*2 + 1][i*2 + 1] = 1
            f[i*2] = f[i*2 + 1] = 0

    return matrix, f


# Решение СЛАУ методом Гаусса
def solve_slae(matrix, f):
    global ULTRACOUNTER

    # Размерность системы
    alen = len(f)
    flag = 0

    # Прямой ход, зануление элементов под диагональю
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
                    #try:
                    f[j] -= m * f[h]
                    #    flag = 0
                    #except:
                    #    flag = 1
                    flag = 0
                except:
                    flag = 1
    result = [0] * alen

    if flag:
        return result, flag

    if -0.001 < matrix[alen - 1][alen - 1] < 0.001:
        return result, flag

    # Обратный ход
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
def update_draft():
    global global_point_list, global_line_list, cPointlist

    cPointlist = []
    cPointlist = global_point_list.copy()
    if update_primitives():
        global_point_list = cPointlist.copy()

    if len(Constraints) != 0:
        print('-'*45)
        print("Наложенные граничения:")
        i = 0
        for con in Constraints:
            print(str(i) + ": " + con.get_description())
            i += 1

    graph_axes.clear()
    graph_axes.set_xlim(XAXES)
    graph_axes.set_ylim(YAXES)

    point_x = []
    point_y = []
    for point in global_point_list:
        point_x.append(point.v_return()[0])
        point_y.append(point.v_return()[1])
    graph_axes.plot(point_x, point_y, 'o', color='black', picker=True, pickradius=5, zorder=2.5)

    for line in global_line_list:
        x, y = line.xy_return()
        graph_axes.plot(x, y, marker='o', picker=True, pickradius=5)

    graph_axes.grid()
    plt.draw()


# Обработка кликов в рабочем поле
def on_click(event):
    """Нажатие по пустой области"""
    global FLAG_CRE, PointCount, global_point_list
    if check_time(datetime.now(), 0):
        if event.button == 3 and FLAG_CRE == 0:  # Создание точки
            add_point(event)
        elif event.button == 3 and FLAG_CRE == 1:  # Создание отрезка вместе с точками
            if PointCount == 0:  # Выбрана первая точка
                add_point(event)
                PointCount = 1
            elif PointCount == 1:  # Выбрана вторая точка
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

        # Координата клика
        mouse_x, mouse_y = event.mouseevent.xdata, event.mouseevent.ydata
        dist = np.sqrt((np.array(xdata) - mouse_x)**2 + (np.array(ydata) - mouse_y)**2)
        mindist = min(dist)

        #distance_to_point = 1.0  # Поиск индекса ближайшей точки
        pointInd = -1
        for point in global_point_list:
            distance_to_point = np.sqrt((point.aCoord[0] - mouse_x)**2 + (point.aCoord[1] - mouse_y)**2)
            if distance_to_point < 0.1:
                pointInd = global_point_list.index(point)
        #distance_to_point = 1.0

        # Если кликнули не по точке
        if mindist > 0.1 and pointInd == -1:
            if FLAG_DEL:  # при удалении линии
                for line in global_line_list:
                    if get_line_by_coord(line, xdata, ydata):
                        delete_line(line)
                        FLAG_DEL = 0
                        button_delete_primitive.color = 'white'

            if FLAG_HOR:  # горизонтальность
                for line in global_line_list:
                    if get_line_by_coord(line, xdata, ydata):
                        tmp = [line.oPoint1, line.oPoint2]
                        tmp2 = constraints.horizontal.Horizontal(tmp)
                        Constraints.append(tmp2)
                        update_draft()
                        FLAG_HOR = 0
                        message_box.set_val('Ограничение наложено')

            if FLAG_VER:  # вертикальность
                for line in global_line_list:
                    if get_line_by_coord(line, xdata, ydata):
                        tmp = [line.oPoint1, line.oPoint2]
                        tmp2 = constraints.vertical.Vertical(tmp)
                        Constraints.append(tmp2)
                        update_draft()
                        FLAG_VER = 0
                        message_box.set_val('Ограничение наложено')

            if FLAG_PAR == 2:  # параллельность
                for line in global_line_list:
                    if get_line_by_coord(line, xdata, ydata):
                        tmp.append(line.oPoint1)
                        tmp.append(line.oPoint2)
                        tmp2 = constraints.parallel.Parallel(tmp)
                        Constraints.append(tmp2)
                        update_draft()
                        FLAG_PAR = 0
                        message_box.set_val('Ограничение наложено')
            if FLAG_PAR == 1:
                for line in global_line_list:
                    if get_line_by_coord(line, xdata, ydata):
                        tmp = [line.oPoint1, line.oPoint2]
                        message_box.set_val('Параллельность. Выберите второй отрезок')
                        FLAG_PAR = 2

            if FLAG_PER == 2:  # перпендикулярность
                for line in global_line_list:
                    if get_line_by_coord(line, xdata, ydata):
                        tmp.append(line.oPoint1)
                        tmp.append(line.oPoint2)
                        tmp2 = constraints.perpendicular.Perpendicular(tmp)
                        Constraints.append(tmp2)
                        update_draft()
                        FLAG_PER = 0
                        message_box.set_val('Ограничение наложено')
            if FLAG_PER == 1:
                for line in global_line_list:
                    if get_line_by_coord(line, xdata, ydata):
                        tmp = [line.oPoint1, line.oPoint2]
                        message_box.set_val('Перпендикулярность. Выберите второй отрезок')
                        FLAG_PER = 2

            if FLAG_POL == 1:  # точка на линии
                for line in global_line_list:
                    if get_line_by_coord(line, xdata, ydata):
                        tmp = [line.oPoint1, line.oPoint2]
                        message_box.set_val('Точка на прямой. Выберите точку')
                        FLAG_POL = 2

            if FLAG_ANG == 2:  # угол, выбор второго отрезка и задание угла
                for line in global_line_list:
                    if get_line_by_coord(line, xdata, ydata):
                        tmp.append(line.oPoint1)
                        tmp.append(line.oPoint2)
                        message_box.set_val('Угол. Введите значение в терминале (градусы)')
                        print('Введите значение угла: ')
                        u = input()
                        flag = 1
                        try:
                            u = float(u)
                        except:
                            flag = 0

                        if flag:
                            tmp2 = constraints.angle.Angle(tmp, u)
                            Constraints.append(tmp2)
                            update_draft()
                            FLAG_ANG = 0
                            message_box.set_val('Ограничение наложено')

                        else:
                            message_box.set_val('Некорректно введено значение угла')
                            print("Некорректно введено значение угла")
                            FLAG_ANG = 0
            if FLAG_ANG == 1:
                for line in global_line_list:
                    if get_line_by_coord(line, xdata, ydata):
                        tmp = [line.oPoint1, line.oPoint2]
                        message_box.set_val('Угол. Выберите второй отрезок')
                        FLAG_ANG = 2

        # Если кликнули на точку
        else:
            points = tuple(zip(xdata[ind], ydata[ind]))
            if FLAG_DEL:  # при удалении точки
                if pointInd != -1:
                    delete_point(global_point_list[pointInd])
                    FLAG_DEL = 0
                    button_delete_primitive.color = 'white'
                    pointInd = -1
                else:
                    for point in global_point_list:
                        if point.v_return() == [round(points[0][0], 10), round(points[0][1], 10)]:
                            delete_point(point)
                            FLAG_DEL = 0
                            button_delete_primitive.color = 'white'

            if FLAG_CRE == 2 and event.mouseevent.button == 3:  # Добавление отрезка по 2-м существующим точкам
                if PointCount == 0:
                    if pointInd != -1:
                        PointInd = global_point_list.index(global_point_list[pointInd])
                        message_box.set_val('Найдена первая точка')
                        PointCount = 1
                        pointInd = -1
                        return 0
                    else:
                        for point in global_point_list:
                            if point.v_return() == [round(points[0][0], 10), round(points[0][1], 10)]:
                                PointInd = global_point_list.index(point)
                                message_box.set_val('Найдена первая точка')
                                PointCount = 1
                                return 0

                if PointCount == 1:
                    if pointInd != -1:
                        if global_point_list[PointInd] != global_point_list[pointInd]:
                            message_box.set_val('Найдена вторая точка')
                            add_line(global_point_list[PointInd], global_point_list[pointInd])
                            PointCount = 0
                            PointInd = 0
                            pointInd = -1
                            return 0
                    else:
                        for point in global_point_list:
                            if point.v_return() == [round(points[0][0], 10), round(points[0][1], 10)]:
                                if global_point_list[PointInd] != point:
                                    message_box.set_val('Найдена вторая точка')
                                    add_line(global_point_list[PointInd], point)
                                    PointCount = 0
                                    PointInd = 0
                                    return 0

            if FLAG_FIX == 1:
                if pointInd != -1:
                    tmp = [global_point_list[pointInd]]
                    tmp2 = constraints.fix.Fixed(tmp, Fixedlist)
                    Constraints.append(tmp2)
                    update_draft()
                    message_box.set_val('Точка зафиксирована')
                    FLAG_FIX = 0
                    pointInd = -1
                else:
                    for point in global_point_list:
                        if point.v_return() == [round(points[0][0], 10), round(points[0][1], 10)]:
                            tmp = [point]
                            tmp2 = constraints.fix.Fixed(tmp, Fixedlist)
                            Constraints.append(tmp2)
                            update_draft()
                            message_box.set_val('Точка зафиксирована')
                            FLAG_FIX = 0
            
            if FLAG_DIS == 2:
                if pointInd != -1:
                    tmp.append(global_point_list[pointInd])
                    message_box.set_val('Расстояние. Введите расстояние в терминале')
                    print('Введите расстояние: ')
                    d = input()
                    flag = 1
                    try:
                        d = float(d)
                    except:
                        flag = 0

                    if flag:
                        tmp2 = constraints.distance.Distance(tmp, d)
                        Constraints.append(tmp2)
                        update_draft()
                        pointInd = -1
                        FLAG_DIS = 0
                        message_box.set_val('Ограничение наложено')
                    else:
                        message_box.set_val('Некорректно введено значение расстояния')
                        print("Некорректно введено значение расстояния")
                        pointInd = -1
                        FLAG_DIS = 0

                else:
                    for point in global_point_list:
                        if point.v_return() == [round(points[0][0], 10), round(points[0][1], 10)]:
                            tmp.append(point)
                            message_box.set_val('Расстояние. Введите расстояние в терминале')
                            print('Введите расстояние: ')
                            d = input()
                            flag = 1

                            try:
                                d = float(d)
                            except:
                                flag = 0

                            if flag:
                                tmp2 = constraints.distance.Distance(tmp, d)
                                Constraints.append(tmp2)
                                update_draft()
                                FLAG_DIS = 0
                            else:
                                message_box.set_val('Некорректно введено значение расстояния')
                                print("Некорректно введено значение расстояния")
                                FLAG_DIS = 0
            if FLAG_DIS == 1:
                if pointInd != -1:
                    tmp = []
                    tmp.append(global_point_list[pointInd])
                    message_box.set_val('Расстояние. Выберите вторую точку')
                    FLAG_DIS = 2
                    pointInd = -1
                else:
                    for point in global_point_list:
                        if point.v_return() == [round(points[0][0], 10), round(points[0][1], 10)]:
                            tmp = []
                            tmp.append(point)
                            message_box.set_val('Расстояние. Выберите вторую точку')
                            FLAG_DIS = 2

            if FLAG_CON == 2:
                if pointInd != -1:
                    tmp.append(global_point_list[pointInd])
                    tmp2 = constraints.coincidence.Coincidence(tmp)
                    Constraints.append(tmp2)
                    update_draft()
                    FLAG_CON = 0
                    pointInd = -1
                    message_box.set_val('Ограничение наложено')
                else:
                    for point in global_point_list:
                        if point.v_return() == [round(points[0][0], 10), round(points[0][1], 10)]:
                            tmp.append(point)
                            tmp2 = constraints.coincidence.Coincidence(tmp)
                            Constraints.append(tmp2)
                            update_draft()
                            FLAG_CON = 0
                            message_box.set_val('Ограничение наложено')

            if FLAG_CON == 1:
                if pointInd != -1:
                    tmp = [global_point_list[pointInd]]
                    message_box.set_val('Совпадение точек. Выберите вторую точку')
                    FLAG_CON = 2
                    pointInd = -1
                else:
                    for point in global_point_list:
                        if point.v_return() == [round(points[0][0], 10), round(points[0][1], 10)]:
                            tmp = [point]
                            message_box.set_val('Совпадение точек. Выберите вторую точку')
                            FLAG_CON = 2

            if FLAG_POL == 2:
                if pointInd != -1:
                    tmp.append(global_point_list[pointInd])
                    tmp2 = constraints.pointOnLine.PointOnLine(tmp)
                    Constraints.append(tmp2)
                    update_draft()
                    FLAG_POL = 0
                    pointInd = -1
                else:
                    for point in global_point_list:
                        if point.v_return() == [round(points[0][0], 10), round(points[0][1], 10)]:
                            tmp.append(point)
                            tmp2 = constraints.pointOnLine.PointOnLine(tmp)
                            Constraints.append(tmp2)
                            update_draft()
                            FLAG_POL = 0


def get_line_by_coord(line, x, y):
    a = line.xy_return()[0][0] == round(x[0], 10)
    b = line.xy_return()[0][1] == round(x[1], 10)
    c = line.xy_return()[1][0] == round(y[0], 10)
    d = line.xy_return()[1][1] == round(y[1], 10)
    return a and b and c and d


# Создание примитивов
def add_point(event):
    # Добавление новой точки
    global global_point_list
    new_point = primitives.point.Point([event.xdata, event.ydata])
    global_point_list.append(new_point)
    update_draft()


def add_line(p1, p2):
    # Добавление нового отрезка
    global global_line_list
    new_line = primitives.line.Line(p1, p2)
    global_line_list.append(new_line)
    update_draft()


# Удаление примитивов
def remove_point_in_constraints(point):
    while True:
        ttt = 0
        for con in Constraints:
            if point in con.Points:
                ttt = 1
                Constraints.remove(con)
        if ttt == 0:
            break

    global_point_list.remove(point)
    point.delPoint()


def delete_point(point_to_delete):
    '''Удаление точки'''
    global global_point_list, global_line_list
    # Если к точке привязаны прямые, то они тоже стираются
    line_list_copy = []
    for line in global_line_list:
        if point_to_delete == line.oPoint1 or point_to_delete == line.oPoint2:
            line_list_copy.append(line)
    for line in line_list_copy:
        global_line_list.remove(line)
        line.delLine()

    remove_point_in_constraints(point_to_delete)

    # Перерисовываем картинку
    update_draft()


def delete_line(line_to_delete):
    '''Удаление отрезка'''
    global global_point_list, global_line_list
    p1 = global_point_list[global_point_list.index(line_to_delete.oPoint1)]
    p2 = global_point_list[global_point_list.index(line_to_delete.oPoint2)]
    p1_incident_count = -1
    p2_incident_count = -1

    # Считаем, со сколькими отрезками смежны точки из удаляемого отрезка
    for oLine_ in global_line_list:
        if p1 == oLine_.oPoint1 or p1 == oLine_.oPoint2:
            p1_incident_count += 1
        if p2 == oLine_.oPoint1 or p2 == oLine_.oPoint2:
            p2_incident_count += 1

    # Удалем прямую и точки
    global_line_list.remove(line_to_delete)
    line_to_delete.delLine()

    # Если точка Линии больше нигде не используется, то её стираем, иначе не трогаем
    if not p1_incident_count:
        remove_point_in_constraints(p1)
    if not p2_incident_count:
        remove_point_in_constraints(p2)
    # Перерисовываем картинку
    update_draft()


# Переключение режимов создания/удаления примитивов
def choose_delete_primitive_mode(event):
    """Обработчик события для кнопки "Удалить точку\""""
    global FLAG_DEL
    if FLAG_DEL == 1:
        FLAG_DEL = 0
        button_delete_primitive.color = 'white'
    else:
        FLAG_DEL = 1
        button_delete_primitive.color = 'grey'


def choose_delete_constraint_mode(event):
    """Обработчик события для кнопки "Удалить ограничение\""""
    global FLAG_DELCON
    if FLAG_DELCON == 1:
        FLAG_DELCON = 0
        button_delete_constraint.color = 'white'
    else:
        if len(Constraints) != 0:
            FLAG_DELCON = 1
            button_delete_constraint.color = 'grey'
            message_box.set_val('Укажите в терминале номер ограничения для удаления')
            print("Номер ограничения для удаления:")
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

            FLAG_DELCON = 0
            button_delete_constraint.color = 'white'
        else:
            print("Ограничений нет, удалять нечего")


def choose_creation_mode(label):
    """Выбор примитива для ввода"""
    global FLAG_CRE
    if radiobuttons_creationtype.value_selected == 'Точка':
        FLAG_CRE = 0
    elif radiobuttons_creationtype.value_selected == 'Отрезок':
        FLAG_CRE = 1
    else:
        FLAG_CRE = 2


# Обработка нажатия на кнопки ограничений
def click_fix_point(event):
    global FLAG_FIX

    if FLAG_FIX == 1:  # Начать заново при повторном клике на кнопку
        FLAG_FIX = 0
    else:
        message_box.set_val("Фиксация. Выберите точку")
        reset_flags()
        FLAG_FIX = 1


def click_distance(event):  # Расстояние
    global FLAG_DIS

    if FLAG_DIS == 1 or FLAG_DIS == 2:
        FLAG_DIS = 0
    else:
        message_box.set_val('Расстояние. Выберите первую точку')
        reset_flags()
        FLAG_DIS = 1


def click_point_on_line(event):  # точка на линии
    global FLAG_POL

    if FLAG_POL == 1 or FLAG_POL == 2:
        FLAG_POL = 0
    else:
        message_box.set_val('Точка на прямой. Выберите отрезок, задающий прямую')
        reset_flags()
        FLAG_POL = 1


def click_coincidence(event):  # совпадение точек
    global FLAG_CON

    if FLAG_CON == 1 or FLAG_CON == 2:
        FLAG_CON = 0
    else:
        message_box.set_val('Совпадение точек. Выберите первую точку')
        reset_flags()
        FLAG_CON = 1


def click_angle(event):  # угол
    global FLAG_ANG

    if FLAG_ANG == 1 or FLAG_ANG == 2:
        FLAG_ANG = 0
    else:
        message_box.set_val('Угол. Выберите первый отрезок')
        reset_flags()
        FLAG_ANG = 1


def click_parallel(event):  # параллельность
    global FLAG_PAR

    if FLAG_PAR == 1 or FLAG_PAR == 2:
        FLAG_PAR = 0
    else:
        message_box.set_val('Параллельность. Выберите первый отрезок')
        reset_flags()
        FLAG_PAR = 1


def click_perpendicular(event):  # перпендикулярность
    global FLAG_PER

    if FLAG_PER == 1 or FLAG_PER == 2:
        FLAG_PER = 0
    else:
        message_box.set_val('Перпендикулярность. Выберите первый отрезок')
        reset_flags()
        FLAG_PER = 1


def click_vertical(event):
    global FLAG_VER

    if FLAG_VER == 1:
        FLAG_VER = 0
    else:
        message_box.set_val('Вертикальность. Выберите отрезок')
        reset_flags()
        FLAG_VER = 1


def click_horizontal(event):
    global FLAG_HOR

    if FLAG_HOR == 1:
        FLAG_HOR = 0
    else:
        message_box.set_val('Горизонтальность. Выберите отрезок')
        reset_flags()
        FLAG_HOR = 1


def init_button(place, title, icon, callback):
    img = None
    label = None
    if icon:
        img = Image.open(icon)
    if title:
        label = title
    button = Button(place, label, img)
    button.on_clicked(callback)
    return button


if __name__ == "__main__":
    fig, graph_axes = plt.subplots()
    graph_axes.set_xlim(XAXES)
    graph_axes.set_ylim(YAXES)
    graph_axes.aspect = 1
    graph_axes.grid()
    plt.tight_layout()

    # Размеры окна и место для кнопок
    fig.set_size_inches(12, 10)
    fig.subplots_adjust(left=0.05, right=0.75, top=0.95, bottom=0.15)

    # Кнопки для ограничений для точек
    button_fix_point = init_button(plt.axes([0.8, 0.865, 0.09, 0.09]), None, "./img/fix.png", click_fix_point)
    button_distance = init_button(plt.axes([0.88, 0.865, 0.09, 0.09]), None, "./img/distance.png", click_distance)
    button_coincidence = init_button(plt.axes([0.8, 0.77, 0.09, 0.09]), None, "./img/coincidence.png", click_coincidence)
    button_point_on_line = init_button(plt.axes([0.88, 0.77, 0.09, 0.09]), None, "./img/point_on_line.png", click_point_on_line)

    # Кнопки ограничений для отрезков
    button_parallel = init_button(plt.axes([0.8, 0.65, 0.09, 0.09]), None, "./img/parallel.png", click_parallel)
    button_perpendicular = init_button(plt.axes([0.88, 0.65, 0.09, 0.09]), None, "./img/perpendicular.png", click_perpendicular)
    button_vertical = init_button(plt.axes([0.8, 0.555, 0.09, 0.09]), None, "./img/vertical.png", click_vertical)
    button_horizontal = init_button(plt.axes([0.88, 0.555, 0.09, 0.09]), None, "./img/horizontal.png", click_horizontal)
    button_angle = init_button(plt.axes([0.8, 0.46, 0.09, 0.09]), None, "./img/angle.png", click_angle)

    # Кнопки для добаления примитивов
    axes_radiobuttons = plt.axes([0.805, 0.31, 0.18, 0.1])
    radiobuttons_creationtype = RadioButtons(axes_radiobuttons, ['Точка', 'Отрезок', 'Отрезок по точкам'])
    radiobuttons_creationtype.on_clicked(choose_creation_mode)

    # Кнопки удаления примитивов или ограничений
    button_delete_primitive = init_button(plt.axes([0.805, 0.22, 0.18, 0.05]), "Удалить примитив", None, choose_delete_primitive_mode)
    button_delete_constraint = init_button(plt.axes([0.805, 0.16, 0.18, 0.05]), "Удалить ограничение", None, choose_delete_constraint_mode)

    # Обработчики кликов
    fig.canvas.mpl_connect('button_press_event', on_click)
    fig.canvas.mpl_connect('pick_event', on_pick)

    # Окно с сообщениями под рабочим полем
    message_box_size = fig.add_axes([0.07, 0.03, 0.85, 0.08])
    message_box = TextBox(message_box_size, None, textalignment="center")

    plt.show()
