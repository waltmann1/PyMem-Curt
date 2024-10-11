#include "montercarlo_vertex.hpp"


int MonteCarloIntegratorVertex::integrate(void)
{
    int accepted_moves = 0;
    for(int vertex_index = 0; vertex_index<_system.Numvertices; vertex_index++)
    {
        //attempt to move the vertex
        double dx, dy, dz;
        if (m_spherical_move)
        {
            double x = _system.vertices[vertex_index].r.x, y =_system.vertices[vertex_index].r.y, z= _system.vertices[vertex_index].r.z;
            double r = sqrt(x * x + y * y + z * z);
            double sigma = sqrt(m_dx * m_dx + m_dy * m_dy + m_dz * m_dz);
            dx = m_rng->gauss_rng(sigma);
            dy = m_rng->gauss_rng(sigma);
            dz = m_rng->gauss_rng(sigma);
            double x_new = x + dx, y_new = y + dy, z_new = z + dz;
            double r_new = sqrt(x_new * x_new + y_new * y_new + z_new * z_new);
            double scale = r / r_new;
            dx = scale * x_new - x;
            dy = scale * y_new - y;
            dz = scale * z_new - z;
        }
        else
        {
            dx = m_dx * (m_rng->drnd() - 0.5);
            dy = m_dy * (m_rng->drnd() - 0.5);
            dz = m_dz * (m_rng->drnd() - 0.5);
        }
        double energy_i = this->ComputeEnergyFromVertex(vertex_index);
        _system.vertices[vertex_index].r.x+= dx;
        _system.vertices[vertex_index].r.y+= dy;
        _system.vertices[vertex_index].r.z+= dz;
        this->update_vertex_normal(vertex_index);
        double energy_f = this->ComputeEnergyFromVertex(vertex_index);
        double delE = energy_f - energy_i;
        accepted_moves++;
        // std::cout << "vertexmove:"<<vertex_index<<std::endl;
        if (!(delE<0.0))
        {
            if(!(m_rng->drnd() < exp(-delE/get_temperature())))
            {
                _system.vertices[vertex_index].r.x-= dx;
                _system.vertices[vertex_index].r.y-= dy;
                _system.vertices[vertex_index].r.z-= dz;
                this->update_vertex_normal(vertex_index);
                accepted_moves--;
                // std::cout<< "regret"<<std::endl;
            }
        }
    }
    return(accepted_moves);
}
