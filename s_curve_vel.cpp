#include <velocity_smooth_plan/s_curve_vel.h>
namespace moying
{
    namespace navigation
    {
        SCurveVel::SCurveVel()
        {
        }
        void SCurveVel::InitParam(int circulation_time_limit, double time_iteration_step)
        {
            circulation_time_limit_ = circulation_time_limit;
            time_iteration_step_ = time_iteration_step;
        }
        void SCurveVel::calSCurve(double s_total, double s_passed, double max_vel, double acc_lim, double jerk, double &vel)
        {
            double t1, t2, t3;
            double v1, v2, v3;
            double s1, s2, s3;
            double a;
            max_vel = fabs(max_vel);
            double fact_max_vel;
            if (max_vel < pow(acc_lim, 2) / jerk)
            {
                v1 = max_vel / 2;
                t1 = sqrt(max_vel / jerk);
                t2 = 0;
            }
            else
            {
                v1 = 1.0 / 2.0 * pow(acc_lim, 2) / jerk;
                t1 = acc_lim / jerk;
                t2 = max_vel / acc_lim - acc_lim / jerk;
            }
            v2 = v1 + acc_lim * t2;
            t3 = t1;
            a = jerk * t1;
            s1 = 1.0 / 6.0 * jerk * pow(t1, 3);
            s2 = 1.0 / 2.0 * a * pow(t2, 2) + v1 * t2;
            s3 = -1.0 / 6.0 * jerk * pow(t3, 3) + 1.0 / 2.0 * a * pow(t3, 2) + v2 * t3;

            if (s_total <= 2 * (s1 + s2 + s3)) //can not reach max_vel
            {
                if (s2 == 0) //acc_lim  Unavailable
                {
                    fact_max_vel = pow((s_total / 2 * pow(jerk, 1.0 / 2.0)), 2.0 / 3.0);
                }
                else
                {
                    if (s_total <= 2 * (jerk * pow(acc_lim / jerk, 3))) //acc_lim available but can not reach
                    {
                        fact_max_vel = pow((s_total / 2 * pow(jerk, 1.0 / 2.0)), 2.0 / 3.0);
                    }
                    else //acc_lim available
                    {
                        double acc_time = (-3 * pow(acc_lim, 2) + pow((pow(acc_lim, 4) + 8 * acc_lim * s_total / 2 * pow(jerk, 2)), 1.0 / 2.0)) / (2 * acc_lim * jerk);
                        fact_max_vel = pow(acc_lim, 2) / jerk + acc_lim * acc_time;
                    }
                }
                if (s_passed < s_total / 2)
                {
                    calStartSCurve(s_passed, fact_max_vel, acc_lim, jerk, vel);
                }
                else if (s_passed >= s_total / 2 && s_passed <= s_total)
                {
                    calBrakeSCurve(s_passed - s_total / 2, fact_max_vel, acc_lim, jerk, vel);
                }
                else
                {
                    vel = 0;
                }
            }
            else //max_vel available
            {
                if (s_passed < (s1 + s2 + s3))
                {
                    calStartSCurve(s_passed, max_vel, acc_lim, jerk, vel);
                }
                else if (s_passed >= (s1 + s2 + s3) && (s_total - s_passed) >= (s1 + s2 + s3))
                {
                    vel = max_vel;
                }
                else if (s_passed > s_total - (s1 + s2 + s3) && s_passed <= s_total)
                {
                    calBrakeSCurve(s_passed + s1 + s2 + s3 - s_total, max_vel, acc_lim, jerk, vel);
                }
                else
                {
                    vel = 0;
                }
            }
        }

        void SCurveVel::calStartSCurve(double s_passed, double end_vel, double acc_lim, double jerk, double &vel)
        {
            double t1, t2, t3;
            double v1, v2, v3;
            double s1, s2, s3;
            double a;
            double max_vel = fabs(end_vel);
            if (max_vel < pow(acc_lim, 2) / jerk)
            {
                v1 = max_vel / 2;
                t1 = sqrt(max_vel / jerk);
                t2 = 0;
            }
            else
            {
                v1 = 1.0 / 2.0 * pow(acc_lim, 2) / jerk;
                t1 = acc_lim / jerk;
                t2 = max_vel / acc_lim - acc_lim / jerk;
            }
            v2 = v1 + acc_lim * t2;
            t3 = t1;
            a = jerk * t1;
            s1 = 1.0 / 6.0 * jerk * pow(t1, 3);
            s2 = 1.0 / 2.0 * a * pow(t2, 2) + v1 * t2;
            s3 = -1.0 / 6.0 * jerk * pow(t3, 3) + 1.0 / 2.0 * a * pow(t3, 2) + v2 * t3;
            if (s_passed <= 0.0007)
            {
                vel = 0.01;
            }
            else if (s_passed < s1 && s_passed > 0.0007)
            {
                vel = 1.0 / 2.0 * jerk * pow(6 * s_passed / jerk, 2.0 / 3.0);
            }
            else if (s_passed >= s1 && s_passed <= (s1 + s2))
            {
                vel = sqrt(pow(v1, 2) + 2 * a * (s_passed - s1));
            }
            else if (s_passed >= (s1 + s2) && s_passed <= (s1 + s2 + s3))
            {
                double t = a / jerk + 0.1;
                double t0, f, fd;
                t0 = 0;
                while (fabs(t - t0) >= 0.001)
                {
                    t0 = t;
                    f = -1.0 / 6.0 * jerk * pow(t0, 3) + 1.0 / 2.0 * a * pow(t0, 2) + v2 * t0 - (s_passed - s1 - s2);
                    fd = -1.0 / 2.0 * jerk * pow(t0, 2) + a * t0 + v2;
                    t = t0 - f / fd;
                }
                vel = -1.0 / 2.0 * jerk * pow(t, 2) + a * t + v2;
            }
            else
            {
                vel = end_vel;
            }
        }

        void SCurveVel::calBrakeSCurve(double s_passed, double begin_vel, double acc_lim, double jerk, double &vel)
        {
            double t1, t2, t3;
            double v1, v2, v3;
            double s1, s2, s3;
            double a;
            double max_vel = fabs(begin_vel);
            if (max_vel < pow(acc_lim, 2) / jerk)
            {
                v1 = max_vel / 2;
                t1 = sqrt(max_vel / jerk);
                t2 = 0;
            }
            else
            {
                v1 = 1.0 / 2.0 * pow(acc_lim, 2) / jerk;
                t1 = acc_lim / jerk;
                t2 = max_vel / acc_lim - acc_lim / jerk;
            }
            v2 = v1 + acc_lim * t2;
            t3 = t1;
            a = jerk * t1;
            s1 = 1.0 / 6.0 * jerk * pow(t1, 3);
            s2 = 1.0 / 2.0 * a * pow(t2, 2) + v1 * t2;
            s3 = -1.0 / 6.0 * jerk * pow(t3, 3) + 1.0 / 2.0 * a * pow(t3, 2) + v2 * t3;

            if (s_passed < s3)
            {
                double t = a / jerk + 0.1;
                double t0, f, fd;
                t0 = 0;
                while (fabs(t - t0) >= 0.001)
                {
                    t0 = t;
                    f = -1.0 / 6.0 * jerk * pow(t0, 3) + 1.0 / 2.0 * a * pow(t0, 2) + v2 * t0 - (s3 - s_passed);
                    fd = -1.0 / 2.0 * jerk * pow(t0, 2) + a * t0 + v2;
                    t = t0 - f / fd;
                }
                vel = -1.0 / 2.0 * jerk * pow(t, 2) + a * t + v2;
            }
            else if (s_passed >= s3 && s_passed <= (s3 + s2))
            {
                vel = sqrt(pow(v1, 2) + 2 * a * (s2 + s3 - s_passed));
            }
            else if (s_passed >= (s3 + s2) && (s1 + s2 + s3) - s_passed > 0.0007)
            {
                vel = 1.0 / 2.0 * jerk * pow(6 * (s1 + s2 + s3 - s_passed) / jerk, 2.0 / 3.0);
            }
            else if ((s1 + s2 + s3) - s_passed <= 0.0007 && (s1 + s2 + s3) - s_passed > 0)
            {
                vel = 0.01;
            }
            else
            {
                vel = 0;
            }
        }

        double SCurveVel::getBrakeDistance(double current_vel, double acc_lim, double jerk, double &t1, double &t2, double &t3)
        {
            double brake_distance = getStartDistance(current_vel, acc_lim, jerk, t1, t2, t3);
            return brake_distance;
        }

        double SCurveVel::getStartDistance(double end_vel, double acc_lim, double jerk, double &t1, double &t2, double &t3)
        {
            double v1, v2, v3;
            double s1, s2, s3;
            double a;
            double start_distance;
            double max_vel = fabs(end_vel);
            if (max_vel < pow(acc_lim, 2) / jerk)
            {
                v1 = max_vel / 2;
                t1 = sqrt(max_vel / jerk);
                t2 = 0;
            }
            else
            {
                v1 = 1.0 / 2.0 * pow(acc_lim, 2) / jerk;
                t1 = acc_lim / jerk;
                t2 = max_vel / acc_lim - acc_lim / jerk;
            }
            v2 = v1 + acc_lim * t2;
            t3 = t1;
            a = jerk * t1;
            s1 = 1.0 / 6.0 * jerk * pow(t1, 3);
            s2 = 1.0 / 2.0 * a * pow(t2, 2) + v1 * t2;
            s3 = -1.0 / 6.0 * jerk * pow(t3, 3) + 1.0 / 2.0 * a * pow(t3, 2) + v2 * t3;
            // cout<<"s1,s2,s3 = "<<s1<<" , "<<s2<<" , "<<s3<<endl;
            start_distance = s1 + s2 + s3;
            return start_distance;
        }

        //scurve_param:[t1 t2 t3, max_vel_time, max_vel]
        bool SCurveVel::getSCurveTime(double s_total, double max_vel, double acc_lim, double jerk, std::vector<double> &scurve_param, double &time)
        {
            if (s_total < 0)
            {
                std::cout << "s_total cannot be negative" << std::endl;
                return false;
            }
            std::vector<double>().swap(scurve_param);

            double t1, t2, t3;
            double v1, v2, v3;
            double s1, s2, s3;
            double a;
            max_vel = fabs(max_vel);
            acc_lim = fabs(acc_lim);
            double fact_max_vel;
            if (max_vel < pow(acc_lim, 2) / jerk)
            {
                v1 = max_vel / 2;
                t1 = sqrt(max_vel / jerk);
                t2 = 0;
            }
            else
            {
                v1 = 1.0 / 2.0 * pow(acc_lim, 2) / jerk;
                t1 = acc_lim / jerk;
                t2 = max_vel / acc_lim - acc_lim / jerk;
            }
            v2 = v1 + acc_lim * t2;
            t3 = t1;
            a = jerk * t1;
            s1 = 1.0 / 6.0 * jerk * pow(t1, 3);
            s2 = 1.0 / 2.0 * a * pow(t2, 2) + v1 * t2;
            s3 = -1.0 / 6.0 * jerk * pow(t3, 3) + 1.0 / 2.0 * a * pow(t3, 2) + v2 * t3;
            double start_dis, start_time, brake_time, max_vel_time;
            if (s_total <= 2 * (s1 + s2 + s3)) //can not reach max_vel
            {
                if (s2 == 0) //acc_lim  Unavailable
                {
                    fact_max_vel = pow((s_total / 2 * pow(jerk, 1.0 / 2.0)), 2.0 / 3.0);
                }
                else
                {
                    if (s_total <= 2 * (jerk * pow(acc_lim / jerk, 3))) //acc_lim available but can not reach
                    {
                        fact_max_vel = pow((s_total / 2 * pow(jerk, 1.0 / 2.0)), 2.0 / 3.0);
                    }
                    else //acc_lim available
                    {
                        double acc_time = (-3 * pow(acc_lim, 2) + pow((pow(acc_lim, 4) + 8 * acc_lim * s_total / 2 * pow(jerk, 2)), 1.0 / 2.0)) / (2 * acc_lim * jerk);
                        fact_max_vel = pow(acc_lim, 2) / jerk + acc_lim * acc_time;
                    }
                }
                start_dis = getStartDistance(fact_max_vel, acc_lim, jerk, t1, t2, t3);
                start_time = t1 + t2 + t3;
                brake_time = start_time;
                max_vel_time = 0;
            }
            else //max_vel available
            {
                fact_max_vel = max_vel;
                start_dis = getStartDistance(fact_max_vel, acc_lim, jerk, t1, t2, t3);
                start_time = t1 + t2 + t3;
                brake_time = start_time;
                max_vel_time = (s_total - 2 * start_dis) / max_vel;
            }
            scurve_param.push_back(t1);
            scurve_param.push_back(t2);
            scurve_param.push_back(t3);
            scurve_param.push_back(max_vel_time);
            scurve_param.push_back(fact_max_vel);
            time = start_time + brake_time + max_vel_time;
            return true;
        }

        //scurve_param:[t1, t2, t3, max_vel_time, max_vel]
        void SCurveVel::calSCurve(double time_passed, double jerk, const vector<double> &scurve_param, double &s_passed, double &vel)
        {
            double t1, t2, t3, max_vel_time, max_vel;
            t1 = scurve_param[0];
            t2 = scurve_param[1];
            t3 = scurve_param[2];
            max_vel_time = scurve_param[3];
            max_vel = scurve_param[4];

            double v1, v2, v3;
            double s1, s2, s3, s4, s_total;
            double a, t_total;
            v1 = 0.5 * jerk * pow(t1, 2);
            a = jerk * t1;
            v2 = v1 + a * t2;
            v3 = max_vel;
            s1 = 1.0 / 6.0 * jerk * pow(t1, 3);
            s2 = 1.0 / 2.0 * a * pow(t2, 2) + v1 * t2;
            s3 = -1.0 / 6.0 * jerk * pow(t3, 3) + 1.0 / 2.0 * a * pow(t3, 2) + v2 * t3;
            s4 = max_vel * max_vel_time;
            s_total = s4 + 2 * (s1 + s2 + s3);
            t_total = max_vel_time + 2 * (t1 + t2 + t3);
            double t;
            if (time_passed > 0 && time_passed <= t1)
            {
                //s = (1/6)jt^3
                t = time_passed;
                s_passed = 1.0 / 6.0 * jerk * pow(t, 3);
                vel = 0.5 * jerk * pow(t, 2);
            }
            else if (time_passed > t1 && time_passed <= t1 + t2)
            {
                //s = s1 + v1*t + 0.5*a*t^2
                t = time_passed - t1;
                s_passed = s1 + v1 * t + 0.5 * a * pow(t, 2);
                vel = v1 + a * t;
            }
            else if (time_passed > t1 + t2 && time_passed <= t1 + t2 + t3)
            {
                //s = s1 + s2 + v2*t + 0.5*a*t^2 - 1/6 *j*t^3
                t = time_passed - t1 - t2;
                s_passed = s1 + s2 + v2 * t + 0.5 * a * pow(t, 2) - 1.0 / 6.0 * jerk * pow(t, 3);
                vel = v2 + a * t - 0.5 * jerk * pow(t, 2);
            }
            else if (time_passed > t1 + t2 + t3 && time_passed <= t1 + t2 + t3 + max_vel_time)
            {
                //s = s1 + s2 + s3 + max_vel*t
                t = time_passed - t1 - t2 - t3;
                s_passed = s1 + s2 + s3 + max_vel * t;
                vel = max_vel;
            }
            else if (time_passed > t1 + t2 + t3 + max_vel_time && time_passed <= t1 + t2 + t3 + max_vel_time + t3)
            {
                //s = s_total - s(t_toal - time_passed)
                t = t_total - time_passed - t1 - t2;
                s_passed = s_total - (s1 + s2 + v2 * t + 0.5 * a * pow(t, 2) - 1.0 / 6.0 * jerk * pow(t, 3));
                vel = v2 + a * t - 0.5 * jerk * pow(t, 2);
            }
            else if (time_passed > t1 + t2 + t3 + max_vel_time + t3 && time_passed <= t1 + t2 + t3 + max_vel_time + t3 + t2)
            {
                //s = s_total - s(t_toal - time_passed)
                t = t_total - time_passed - t1;
                s_passed = s_total - (s1 + v1 * t + 0.5 * a * pow(t, 2));
                vel = v1 + a * t;
            }
            else if (time_passed > t1 + t2 + t3 + max_vel_time + t3 + t2 && time_passed <= t1 + t2 + t3 + max_vel_time + t3 + t2 + t1)
            {
                //s = s_total - s(t_toal - time_passed)
                t = t_total - time_passed;
                s_passed = s_total - (1.0 / 6.0 * jerk * pow(t, 3));
                vel = 0.5 * jerk * pow(t, 2);
            }
            else if (time_passed > t_total)
            {
                //s = s_total
                s_passed = s_total;
                vel = 0;
            }
        }

        bool SCurveVel::getMaxVel(double s_total, double time_total, double jerk, double acc_lim, double &max_vel, vector<double> &scurve_param)
        {
            double time;
            if (max_vel <= 0)
                max_vel = 0.1;
            if (s_total < 0.000001)
            {
                max_vel = 0;
                std::vector<double>().swap(scurve_param);
                for (int i = 0; i < 5; i++)
                    scurve_param.push_back(0.0);
                return true;
            }

            if (!getSCurveTime(s_total, max_vel, acc_lim, jerk, scurve_param, time))
                return false;
            int64_t time_begin = GetTimeStamp();
            int64_t time_now;
            double up_boundary, down_boundary, test_value;

            if (time_total == time)
            {
                return true;
            }
            else if (time_total > time) //The vel should be decreased
            {
                up_boundary = scurve_param[4];
                down_boundary = scurve_param[4];
                while (1 > 0)
                {
                    down_boundary = down_boundary / 2.0;
                    time_now = GetTimeStamp();
                    if (time_now - time_begin > circulation_time_limit_)
                    {
                        cout << "SCurveVel::getMaxVel:   Calculation takes too long to find down boundary" << endl;
                        cout << "time_now - time_begin = " << time_now - time_begin << endl;
                        return false;
                    }
                    if (!getSCurveTime(s_total, down_boundary, acc_lim, jerk, scurve_param, time))
                        return false;
                    if (time > time_total)
                        break;
                }
            }
            else //The vel should be increased
            {
                down_boundary = scurve_param[4];
                up_boundary = scurve_param[4];
                while (1 > 0)
                {
                    up_boundary = 2 * up_boundary;
                    time_now = GetTimeStamp();
                    if (time_now - time_begin > circulation_time_limit_)
                    {
                        cout << "SCurveVel::getMaxVel:   Calculation takes too long to find up boundary" << endl;
                        cout << "time_now - time_begin = " << time_now - time_begin << endl;
                        return false;
                    }
                    if (!getSCurveTime(s_total, up_boundary, acc_lim, jerk, scurve_param, time))
                        return false;
                    if (time < time_total)
                        break;
                }
            }
            // cout << "boundary: (up , down) = (" << up_boundary << " , " << down_boundary << ")" << endl;
            while (1 > 0)
            {
                time_now = GetTimeStamp();
                if (time_now - time_begin > circulation_time_limit_)
                {
                    cout << "SCurveVel::getMaxVel:   Calculation takes too long to find time" << endl;
                    cout << "time_now - time_begin = " << time_now - time_begin << endl;
                    return false;
                }
                test_value = (up_boundary + down_boundary) / 2.0;
                if (!getSCurveTime(s_total, test_value, acc_lim, jerk, scurve_param, time))
                    return false;
                if (fabs(time - time_total) < time_iteration_step_)
                    break;
                else if (time > time_total)
                    down_boundary = test_value;
                else
                    up_boundary = test_value;
                // cout<<"boundary: (up , down) = ("<<up_boundary<<" , "<<down_boundary<<")"<<endl;
            }
            // cout << "time cost = " << (time_now - time_begin) << endl;
            max_vel = test_value;
            return true;
        }

        int64_t SCurveVel::GetTimeStamp()
        {
            std::chrono::time_point<std::chrono::system_clock, std::chrono::milliseconds> tp = std::chrono::time_point_cast<std::chrono::milliseconds>(std::chrono::system_clock::now());
            auto tmp = std::chrono::duration_cast<std::chrono::milliseconds>(tp.time_since_epoch());
            int64_t timestamp = tmp.count();
            return timestamp;
        }
    }
}
