#ifndef VELOCITY_SMOOTH_PLAN_H
#define VELOCITY_SMOOTH_PLAN_H

#include <math.h>
#include <iostream>
#include <vector>
#include <chrono>
using namespace std;
namespace moying
{
    namespace navigation
    {
        /**
        * @brief A class to calculate  track follow velocity
        */
        class SCurveVel
        {
        public:
            SCurveVel();
            /**
            * @brief param init function
            */
            void InitParam(int circulation_time_limit, double time_iteration_step);
            /**
            * @brief Get adaptive velocity
            * @param[in] s_total , velocity planning total distance
            * @param[in] s_passed , distance traveled
            * @param[in] max_vel , maximum velocity limit
            * @param[in] acc_lim , maximum acceleration limit
            * @param[in] jerk ,jerk  value
            * @param[out] vel, velocity
            */
            void calSCurve(double s_total, double s_passed, double max_vel, double acc_lim, double jerk, double &vel);

            /**
            * @brief Get start velocity
            * @param[in] s_passed , distance traveled
            * @param[in] end_vel , velocity after start
            * @param[in] acc_lim , maximum acceleration limit
            * @param[in] jerk ,jerk  value
            * @param[out] vel, velocity
            */
            void calStartSCurve(double s_passed, double end_vel, double acc_lim, double jerk, double &vel);

            /**
            * @brief Get brake velocity
            * @param[in] s_passed , distance traveled
            * @param[in] begin_vel , velocity before brake
            * @param[in] acc_lim , maximum acceleration limit
            * @param[in] jerk ,jerk  value
            * @param[out] vel, velocity
            */
            void calBrakeSCurve(double s_passed, double begin_vel, double acc_lim, double jerk, double &vel);

            /**
            * @brief Get brake distance 
            * @param[in] current_vel , velocity before brake
            * @param[in] acc_lim , maximum acceleration limit
            * @param[in] jerk ,jerk  value
            * @return brake distance
            */
            double getBrakeDistance(double current_vel, double acc_lim, double jerk, double &t1, double &t2, double &t3);

            /**
            * @brief Get start distance 
            * @param[in] end_vel , velocity after start
            * @param[in] acc_lim , maximum acceleration limit
            * @param[in] jerk ,jerk  value
            * @return start distance
            */
            double getStartDistance(double end_vel, double acc_lim, double jerk,double &t1, double &t2, double &t3);

            bool getSCurveTime(double s_total, double max_vel, double acc_lim, double jerk,std::vector<double> &scurve_param, double &time);
            void calSCurve(double time_passed, double jerk, const vector<double> &scurve_param, double &s_passed, double &vel);
            bool getMaxVel(double s_total, double time_total,double jerk,double acc_lim, double &max_vel,vector<double> &scurve_param);
            int64_t GetTimeStamp();

        private:
            int circulation_time_limit_;
            double time_iteration_step_;
        };
    }
}
#endif