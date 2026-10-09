#pragma once

#include "Bodies.h"

namespace motus {
    
struct GravityParams {
    double G  = 1.0;
    double eps = 1e-3;
};

void computeGravity(const BodiesView& b, const GravityParams& p, Accel& out);
double totalEnergy(const BodiesView& b, const GravityParams& p);


class Gravity {
    public:
        //Constructors
        Gravity() = default;
        Gravity(double G, double eps);

        Gravity(const Gravity&) = default;
        Gravity& operator=(const Gravity&) = default;

        Gravity(Gravity&&) = default;
        Gravity& operator=(Gravity&&) = default;

        //Destructors
        ~Gravity() = default;

        //Core functions
        void compute(BodiesView b, AccelView out) const;
        [[nodiscard]] double potentialEnergy(BodiesView b) const;
        
        //Setters
        void setG       (double newG);
        void resetG     ();
        void setEps     (double newEps);
        void resetEps   ();

        //Getters
        [[nodiscard]] double G()    const noexcept;
        [[nodiscard]] double eps()  const noexcept;


    private:
        static constexpr double     s_defaultG      = 6.674e-11;
        static constexpr double     s_defaultEps    = 1e-3;
        double                      m_G             = s_defaultG;
        double                      m_eps           = s_defaultEps;
};

}
