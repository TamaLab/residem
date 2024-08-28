# script adapted from https://www.geeksforgeeks.org/how-to-check-if-a-given-point-lies-inside-a-polygon/

class PolyPoint:
    """Finds if a point is inside a polygon, the required inputs are [XYZ] and an arry of points [[X1Y1Z1],
    [X2Y2Z3][...]] """
    def __init__(self, point, polygon) -> bool:
        self.point = point
        self.polygon = polygon
        self.int_max = 1000
        self.inside_polygon = self.is_inside_polygon()

    @staticmethod
    def on_segment(p, q, r):
        if ((q[0] <= max(p[0], r[0])) &
                (q[0] >= min(p[0], r[0])) &
                (q[1] <= max(p[1], r[1])) &
                (q[1] >= min(p[1], r[1]))):
            return True
        return False

    @staticmethod
    def orientation(p: tuple, q: tuple, r: tuple) -> int:
        val = (((q[1] - p[1]) *
                (r[0] - q[0])) -
               ((q[0] - p[0]) *
                (r[1] - q[1])))
        if val == 0:
            return 0
        if val > 0:
            return 1
        else:
            return 2


    def do_intersect(self,p1, q1, p2, q2):
        o1 = self.orientation(p1, q1, p2)
        o2 = self.orientation(p1, q1, q2)
        o3 = self.orientation(p2, q2, p1)
        o4 = self.orientation(p2, q2, q1)
        if (o1 != o2) and (o3 != o4):
            return True
        if (o1 == 0) and (self.on_segment(p1, p2, q1)):
            return True
        if (o2 == 0) and (self.on_segment(p1, q2, q1)):
            return True
        if (o3 == 0) and (self.on_segment(p2, p1, q2)):
            return True
        if (o4 == 0) and (self.on_segment(p2, q1, q2)):
            return True
        return False

    def is_inside_polygon(self):
        n = len(self.polygon)
        if n < 3:
            return False
        extreme = (self.int_max, self.point[1])
        decrease = 0
        count = i = 0
        while True:
            next = (i + 1) % n
            if self.polygon[i][1] == self.point[1]:
                decrease += 1
            if (self.do_intersect(self.polygon[i],
                                  self.polygon[next],
                                  self.point, extreme)):
                if self.orientation(self.polygon[i], self.point,
                                    self.polygon[next]) == 0:
                    return self.on_segment(self.polygon[i], self.point,
                                           self.polygon[next])
                count += 1
            i = next
            if i == 0:
                break
        count -= decrease
        return count % 2 == 1




