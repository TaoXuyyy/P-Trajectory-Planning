#!/usr/bin/python3
import numpy as np
from matplotlib import pyplot as plt

class ScurvePlanner():

    def __init__(self,q0,q1,v0,v1,v_max,a_max,j_max):
        q0 = q0
        q1 = q1
        v0 = v0
        v1 = v1
        v_max = v_max
        a_max = a_max
        j_max = j_max

    def _compute_maximum_speed_reached_or_not_reached(self,q0, q1, v0, v1, v_max, a_max, j_max):
        if (v_max-v0)*j_max < a_max**2:               # a_max is not reached
            if v0 > v_max:
                Tj1 = 0
                Ta = 0
                alima = 0

            else:
                Tj1 = np.sqrt((v_max-v0)/j_max)
                Ta = 2*Tj1
                alima = Tj1 * j_max

        else:
            # a_max is reached
            Tj1 = a_max/j_max
            Ta = Tj1 + (v_max-v0)/a_max
            alima = a_max

        # Deceleration period
        if (v_max-v1)*j_max < a_max**2:
            # a_min is not reached
            Tj2 = np.sqrt((v_max-v1)/j_max)
            Td = 2*Tj2
            alimd = Tj2 * j_max

        else:
            # a_min is reached
            Tj2 = a_max/j_max
            Td = Tj2 + (v_max-v1)/a_max
            alimd = a_max
        # print(Ta,Td)
            # 得出Ta和Td，下面要对这两个时间进行判定

        Tv = (q1-q0)/v_max - (Ta/2)*(1+v0/v_max)-(Td/2)*(1+v1/v_max)
        #计算匀速段时间
        # print (Tv)
        T = Tv + Ta + Td

        if Tv > 0:
            vlim = v_max
            T = Tv + Ta + Td
        else:
            Tv = 0
            amax_org = a_max
            delta = (a_max**4)/(j_max**2) + 2*(v0**2+v1**2) + a_max*(4*(q1 - q0) - 2*a_max/j_max*(v0 + v1))
            Tj1 = a_max/j_max
            Ta = (a_max**2/j_max - 2*v0 + np.sqrt(delta))/2/a_max
            Tj2 = a_max/j_max
            Td = ((a_max**2)/j_max - 2*v1 + np.sqrt(delta))/2/a_max
            vlim = v0 + (Ta - Tj1)*alima
            count = 0
            while Ta < 2*Tj1 or Td < 2*Tj2:
                count += 1
                a_max = a_max - amax_org*0.1
                alima = a_max
                alimd = a_max
                if a_max > 0:
                    delta = (a_max ** 4) / (j_max ** 2) + 2 * (v0 ** 2 + v1 ** 2) + a_max * (4 * (q1 - q0) - 2 * a_max / j_max * (v0 + v1))
                else:
                    delta = (a_max ** 4) / (j_max ** 2) + 2 * (v0 ** 2 + v1 ** 2) - a_max * (4 * (q1 - q0) - 2 * a_max / j_max * (v0 + v1))

                Tj1 = a_max/j_max
                Ta = ((a_max**2)/j_max - 2*v0 + np.sqrt(delta))/2/a_max
                Tj2 = a_max/j_max
                Td = ((a_max**2)/j_max - 2*v1 + np.sqrt(delta))/2/a_max
                vlim = v0 + (Ta - Tj1)*alima
                vlima = vlim
                vlimb = v1 - (Td - Tj2)*alimd
            # print (Tj1,Ta,Td,a_max)

            if Ta < 0 or Td < 0:
                if v0 > v1:
                    Ta = 0
                    Tj1 = 0
                    alima = 0
                    Td = 2*(q1 - q0)/(v1 + v0)
                    Tj2 = (j_max*(q1 - q0) - np.sqrt(j_max*(j_max*((q1 - q0)**2) + ((v1 + v0)**2)*(v1 - v0))))/j_max/(v1 + v0)
                    alimd = -j_max*Tj2
                    vlim = v1 - (Td - Tj2)*alimd
                    alimd = -alimd
                else:
                    Td = 0
                    Tj2 = 0

                    Ta = 2*(q1 - q0)/(v1 + v0)
                    Tj1  = (j_max*(q1 - q0) - np.sqrt(j_max*(j_max*((q1 - q0)**2) - ((v1 + v0)**2)*(v1 - v0))))/j_max/(v1 +v0)
                    alima = j_max*Tj1
                    vlim = v0 + (Ta - Tj1)*alima

        self.T = Tv + Ta + Td
        self.Tj1 = Tj1
        self.Tj2 = Tj2
        self.Ta = Ta
        self.Tv = Tv
        self.Td = Td
        self.alima = alima
        self.alimd = alimd
        self.vlim = vlim
        return self.T

    def _get_trajectory_func(self,q0,q1,v0,v1,v_max,a_max,j_max):

        step = 0.001
        t_list = np.arange(0, self.T, step)

        q_list = []
        v_list = []
        a_list = []
        j_list = []

        for t in t_list:
            if 0 <= t < self.Tj1:     #0<= t <0.06
                q = q0 + v0*t + j_max*(t**3)/6
                qd = v0 + j_max * (t ** 2) / 2
                qdd = j_max * t
                qddd = j_max

            elif self.Tj1 <= t < (self.Ta - self.Tj1):   #0.06<= t <0.066
                q = q0 + v0*t + self.alima/6*(3*(t**2) - 3*self.Tj1*t + self.Tj1**2)
                qd = v0 + self.alima * (t - self.Tj1 / 2)
                qdd = self.alima
                qddd = 0

            elif (self.Ta-self.Tj1) <= t < self.Ta:   #0.066<= t <0.126
                q = q0 + (self.vlim+v0)*self.Ta/2 - self.vlim*(self.Ta - t) + j_max*(self.Ta - t)**3/6
                qd = self.vlim - j_max * ((self.Ta - t) ** 2) / 2
                qdd = j_max * (self.Ta - t)
                qddd = -j_max


            # Constant velocity phase
            elif self.Ta <= t < (self.Ta + self.Tv):   #0.126 <= t < 0.126
                q = q0 + (self.vlim+v0)*self.Ta/2 + self.vlim*(t-self.Ta)
                qd = self.vlim
                qdd = 0
                qddd = 0

            # Deceleration phase
            elif (self.T - self.Td) <= t < (self.T-self.Td+self.Tj2):   #0.126 <= t <0.186
                q = q1 - (self.vlim+v1)*self.Td/2 + self.vlim*(t-self.T+self.Td) -j_max*((t-self.T+self.Td)**3)/6
                qd = self.vlim - j_max * ((t-self.T+self.Td) ** 2) / 2
                qdd = -j_max * (t-self.T+self.Td)
                qddd = -j_max


            elif (self.T-self.Td+self.Tj2) <= t < (self.T-self.Tj2): #0.186 <= t < 0.192

                q = q1 - (self.vlim+v1)*self.Td/2 + self.vlim*(t-self.T+self.Td) -self.alimd/6*(3*((t-self.T+self.Td)**2) - 3*self.Tj2*(t-self.T+self.Td) + self.Tj2**2)
                qd = self.vlim - self.alimd * ((t-self.T+self.Td) - self.Tj2 / 2)
                qdd = -self.alimd
                qddd = 0

            elif (self.T-self.Tj2) <= t <= self.T:#0.192 <= t <= 0.252
                q = q1 - v1*(self.T-t) - j_max*((self.T-t)**3)/6
                qd = v1 + j_max * ((self.T-t) ** 2) / 2
                qdd = -j_max * (self.T-t)
                qddd = j_max
            # else:
            #     qdd = 0
            #     qd= v1
            #     q = q1
            #     qddd = 0


            q_list.append(q)
            v_list.append(qd)
            a_list.append(qdd)
            j_list.append(qddd)

        # with open("data.csv", 'w+') as f:
        #     f.write(str(t_list))
        #     f.close()
        return q_list, v_list, a_list, j_list, t_list

if __name__ == "__main__":
    q0 = 0
    q1 = 0.1
    v0 = 0
    v1 = 0
    v_max = 0.1
    a_max = 1
    j_max = 10
    p = ScurvePlanner(q0,q1,v0,v1,v_max,a_max,j_max)
    p._compute_maximum_speed_reached_or_not_reached(q0, q1, v0, v1, v_max, a_max, j_max)
    c = p._get_trajectory_func(q0,q1,v0,v1,v_max,a_max,j_max)

    plt.suptitle("s-curve")
    #
    # # 绘制位置曲线
    plt.subplot(411)
    plt.plot(c[4], c[0], "r")
    plt.ylabel("Position")
    plt.grid('on')
    #
    #
    # # 绘制速度曲线
    plt.subplot(412)
    plt.plot(c[4], c[1], "g")
    plt.ylabel("Velocity")
    plt.grid('on')
    #
    #
    # # 绘制加速度曲线
    plt.subplot(413)
    plt.plot(c[4], c[2], "b")
    plt.ylabel("Acceleration")
    plt.grid('on')
    #
    #
    plt.subplot(414)
    plt.plot(c[4], c[3], "b")
    plt.ylabel("jerk")
    plt.grid('on')

    plt.xlabel("time")
    plt.show()















