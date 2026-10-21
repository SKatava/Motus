#pragma once

#include <memory>

#include "Force.h"
#include "Integrator.h"
#include "Profiler.h"

namespace motus {

class Simulation {
    public:
        //Constructors
        explicit Simulation(IntegratorMethod m = IntegratorMethod::VERLET);

        Simulation(const Simulation&) = delete;
        Simulation& operator=(const Simulation&) = delete;

        Simulation(Simulation&&) = default;
        Simulation& operator=(Simulation&&) = default;

        //Destructors
        ~Simulation() = default;

        //Core functions
        void step(BodiesView, double dt);
        double totalEnergy(BodiesView b);
        double kineticEnergy(BodiesView b);

        //Setters
        
        
        //Getters
        [[nodiscard]]       Gravity& gravity() noexcept;
        [[nodiscard]] const Gravity& gravity() const noexcept;

    private:
        ForceSet                    m_forceSet;
        std::unique_ptr<Integrator> m_integrator; 
        double                      m_time          = 0.0;

};

}
