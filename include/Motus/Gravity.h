#pragma once

#include "Bodies.h"

namespace motus {

class Gravity {
    public:
        enum class Algorithm {
            DIRCET_NEWTON,
            BARNES_HUT
        };

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
        void compute(BodiesView b, std::span<glm::dvec3> out) const;
        [[nodiscard]] double potentialEnergy(BodiesView b) const;
        
        //Setters
        void setG           (double newG);
        void resetG         ();
        void setEps         (double newEps);
        void resetEps       ();
        void setAlgorithm   (Algorithm newAlgo);
        void setTheta       (double newTheta);
        void resetTheta     ();

        //Getters
        [[nodiscard]] double G()        const noexcept;
        [[nodiscard]] double eps()      const noexcept;
        [[nodiscard]] double theta()    const noexcept;

    private:
        static constexpr double     s_defaultG      = 6.674e-11;
        static constexpr double     s_defaultEps    = 1e-3;
        static constexpr double     s_defaultTheta  = 0.5;
        double                      m_G             = s_defaultG;
        double                      m_eps           = s_defaultEps;
        double                      m_theta         = s_defaultTheta;
        Algorithm                   m_algo          = Algorithm::DIRCET_NEWTON;
        

        //Private function
        void directNewton(BodiesView b, AccelView out) const;
        void barnesHut(BodiesView b, AccelView out) const;
};

}
