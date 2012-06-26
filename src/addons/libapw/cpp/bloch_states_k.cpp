#include "lapw.h"

bloch_states_k::bloch_states_k(int ngk) : ngk(ngk)
{
    idxg.resize(ngk);
    idxgfft.resize(ngk);

    lapw_basis_size = ngk + lapw_global.size_wfmt_lo;

    wave_function_size = ngk + lapw_global.size_wfmt;
}

/// copy matching coefficients
void bloch_states_k::copy_apwalm(complex16 *apwalm_)
{   
    apwalm.set_dimensions(ngk, lapw_global.size_wfmt_apw);
    apwalm.allocate();
    
    mdarray<complex16,4> apwalm_tmp(apwalm_, lapw_global.ngkmax, lapw_global.apwordmax, lapw_global.lmmaxapw, lapw_global.atoms.size());
    
    for (int ias = 0; ias < (int)lapw_global.atoms.size(); ias++)
    {
        Atom *atom = lapw_global.atoms[ias];
        Species *species = atom->species;
        
        for (int j = 0; j < species->index.apw_size(); j++)
        {
            int io = species->index[j].order;
            int lm = species->index[j].lm;
            for (int ig = 0; ig < ngk; ig++)
                apwalm(ig, atom->offset_apw + j) = conj(apwalm_tmp(ig, io, lm, ias));
        }
    }
}

inline void move_apw_blocks(complex16 *wf)
{
    for (int ias = (int)lapw_global.atoms.size() - 1; ias > 0; ias--)
    {
        int final_block_offset = lapw_global.atoms[ias]->offset_wfmt;
        int initial_block_offset = lapw_global.atoms[ias]->offset_apw;
        int block_size = lapw_global.atoms[ias]->species->index.apw_size();

        memmove(&wf[final_block_offset], &wf[initial_block_offset], block_size * sizeof(complex16));
    }
}

inline void copy_lo_blocks(complex16 *wf, complex16 *evec)
{
    for (int ias = 0; ias < (int)lapw_global.atoms.size(); ias++)
    {
        int final_block_offset = lapw_global.atoms[ias]->offset_wfmt + lapw_global.atoms[ias]->species->index.apw_size();
        int initial_block_offset = lapw_global.atoms[ias]->offset_lo;
        int block_size = lapw_global.atoms[ias]->species->index.lo_size();
        
        if (block_size > 0)
            memcpy(&wf[final_block_offset], &evec[initial_block_offset], block_size * sizeof(complex16));
    }
}

inline void copy_pw_block(int ngk, complex16 *wf, complex16 *evec)
{
    memcpy(wf, evec, ngk * sizeof(complex16));
}

void bloch_states_k::generate_scalar_wave_functions()
{
    timer t("bloch_states_k::generate_scalar_wave_functions");
    
    scalar_wave_functions.set_dimensions(wave_function_size, lapw_global.nstfv);
    scalar_wave_functions.allocate();
    
    zgemm<cpu>(2, 0, lapw_global.size_wfmt_apw, lapw_global.nstfv, ngk, zone, &apwalm(0, 0), ngk, &evecfv(0, 0), 
        evecfv.size(0), zzero, &scalar_wave_functions(0, 0), scalar_wave_functions.size(0));
    
    for (int j = 0; j < lapw_global.nstfv; j++)
    {
        move_apw_blocks(&scalar_wave_functions(0, j));

        copy_lo_blocks(&scalar_wave_functions(0, j), &evecfv(ngk, j));

        copy_pw_block(ngk, &scalar_wave_functions(lapw_global.size_wfmt, j), &evecfv(0, j));
    }
}

void bloch_states_k::generate_spinor_wave_functions(diagonalization mode)
{
    timer t("bloch_states_k::generate_spinor_wave_functions");
    
    spinor_wave_functions.set_dimensions(wave_function_size, lapw_global.nspinor, lapw_global.nstsv);
    spinor_wave_functions.allocate();

    if (mode == second_variational)
    {
        for (int ispn = 0; ispn < lapw_global.nspinor; ispn++)
            zgemm<cpu>(0, 0, spinor_wave_functions.size(0), lapw_global.nstsv, lapw_global.nstfv, zone, &scalar_wave_functions(0, 0), 
                       scalar_wave_functions.size(0), &evecsv(ispn * lapw_global.nstfv, 0), evecsv.size(0), zzero, 
                       &spinor_wave_functions(0, ispn, 0), spinor_wave_functions.size(0) * spinor_wave_functions.size(1));
    }

    if (mode == full)
    {
        for (int ispn = 0; ispn < lapw_global.nspinor; ispn ++)
        {
            zgemm<cpu>(2, 0, lapw_global.size_wfmt_apw, lapw_global.nstsv, ngk, zone, &apwalm(0, 0), ngk, &evecfd(ispn * lapw_basis_size, 0), 
                       evecfd.size(0), zzero, &spinor_wave_functions(0, ispn, 0), spinor_wave_functions.size(0) * spinor_wave_functions.size(1));
    
            for (int j = 0; j < lapw_global.nstsv; j++)
            {
                move_apw_blocks(&spinor_wave_functions(0, ispn, j));
                
                copy_lo_blocks(&spinor_wave_functions(0, ispn, j), &evecfd(ispn * lapw_basis_size + ngk, j));

                copy_pw_block(ngk, &spinor_wave_functions(lapw_global.size_wfmt, ispn, j), &evecfd(ispn * lapw_basis_size, j));
            }
        }
    }
}

void bloch_states_k::test_scalar_wave_functions(int use_fft)
{
    std::vector<complex16> zfft;
    std::vector<complex16> v1;
    std::vector<complex16> v2;
    
    if (use_fft == 1) 
    {
        v1.resize(ngk);
        zfft.resize(lapw_global.ngrtot);
    }
    
    if (use_fft == 2) 
    {
        v1.resize(lapw_global.ngrtot);
        v2.resize(lapw_global.ngrtot);
    }
    
    double maxerr = 0;

    for (int j1 = 0; j1 < lapw_global.nstfv; j1++)
    {
        if (use_fft == 1)
        {
            memset(&zfft[0], 0, lapw_global.ngrtot * sizeof(complex16));
            for (int ig = 0; ig < ngk; ig++) 
                zfft[idxgfft[ig]] = scalar_wave_functions(lapw_global.size_wfmt + ig, j1);
                
            lapw_fft(1, &zfft[0]);
            
            for (int ir = 0; ir < lapw_global.ngrtot; ir++)
                zfft[ir] *= lapw_global.cfunir[ir];
            
            lapw_fft(-1, &zfft[0]);
            
            for (int ig = 0; ig < ngk; ig++) 
                v1[ig] = zfft[idxgfft[ig]];
        }
        
        if (use_fft == 2)
        {
            memset(&v1[0], 0, lapw_global.ngrtot * sizeof(complex16));
            for (int ig = 0; ig < ngk; ig++) 
                v1[idxgfft[ig]] = scalar_wave_functions(lapw_global.size_wfmt + ig, j1);
            
            lapw_fft(1, &v1[0]);
        }
       
        for (int j2 = 0; j2 < lapw_global.nstfv; j2++)
        {
            complex16 zsum(0,0);
            for (int ias = 0; ias < (int)lapw_global.atoms.size(); ias++)
            {
                int offset_wfmt = lapw_global.atoms[ias]->offset_wfmt;
                Species *sp = lapw_global.atoms[ias]->species;

                for (int l = 0; l <= lapw_global.lmaxapw; l++)
                {
                    int ordmax = sp->radial_index.nrf(l);
                    for (int io1 = 0; io1 < ordmax; io1++)
                        for (int io2 = 0; io2 < ordmax; io2++)
                            for (int m = -l; m <= (int)l; m++)
                                zsum += conj(scalar_wave_functions(offset_wfmt + sp->index(l, m, io1), j1)) *
                                             scalar_wave_functions(offset_wfmt + sp->index(l, m, io2), j2) * 
                                             lapw_runtime.ovlprad(l, io1, io2, ias);
                }
            }
            
            if (use_fft == 0) 
            {
                for (int ig1 = 0; ig1 < ngk; ig1++)
                {
                    for (int ig2 = 0; ig2 < ngk; ig2++)
                    {
                        int ig3 = idxG12(idxg[ig1], idxg[ig2]);
                        zsum += conj(scalar_wave_functions(lapw_global.size_wfmt + ig1, j1)) * scalar_wave_functions(lapw_global.size_wfmt + ig2, j2) * lapw_global.cfunig[ig3];
                    }
               }
           }
           if (use_fft == 1)
           {
               for (int ig = 0; ig < ngk; ig++)
                   zsum += conj(v1[ig]) * scalar_wave_functions(lapw_global.size_wfmt + ig, j2);
           }
           
           if (use_fft == 2)
           {
               memset(&v2[0], 0, lapw_global.ngrtot * sizeof(complex16));
               for (int ig = 0; ig < ngk; ig++) 
                   v2[idxgfft[ig]] = scalar_wave_functions(lapw_global.size_wfmt + ig, j2);
            
               lapw_fft(1, &v2[0]);

               for (int ir = 0; ir < lapw_global.ngrtot; ir++)
                   zsum += conj(v1[ir]) * v2[ir] * lapw_global.cfunir[ir] / double(lapw_global.ngrtot);
           }

           zsum = (j1 == j2) ? zsum - zone : zsum;
           maxerr = std::max(maxerr, abs(zsum));
        }
    }
    std :: cout << "maximum error = " << maxerr << std::endl;
}

void bloch_states_k::test_spinor_wave_functions(int use_fft)
{
    std::vector<complex16> zfft(lapw_global.ngrtot);
    std::vector<complex16> v1[lapw_global.nspinor];
    
    for (int ispn = 0; ispn < lapw_global.nspinor; ispn++)
    {
        if (use_fft == 1) v1[ispn].resize(ngk);
        if (use_fft == 2) v1[ispn].resize(lapw_global.ngrtot);
    }
    
    double maxerr = 0;

    for (int j1 = 0; j1 < lapw_global.nstsv; j1++)
    {
        if (use_fft == 1)
        {
            for (int ispn = 0; ispn < lapw_global.nspinor; ispn++)
            {
                memset(&zfft[0], 0, lapw_global.ngrtot * sizeof(complex16));
                for (int ig = 0; ig < ngk; ig++) 
                    zfft[idxgfft[ig]] = spinor_wave_functions(lapw_global.size_wfmt + ig, ispn, j1);
                    
                lapw_fft(1, &zfft[0]);
                
                for (int ir = 0; ir < lapw_global.ngrtot; ir++)
                    zfft[ir] *= lapw_global.cfunir[ir];
                
                lapw_fft(-1, &zfft[0]);
                
                for (int ig = 0; ig < ngk; ig++) 
                    v1[ispn][ig] = zfft[idxgfft[ig]];
            }
        }
        
        if (use_fft == 2)
        {
            for (int ispn = 0; ispn < lapw_global.nspinor; ispn++)
            {
                memset(&v1[ispn][0], 0, lapw_global.ngrtot * sizeof(complex16));
                for (int ig = 0; ig < ngk; ig++) 
                    v1[ispn][idxgfft[ig]] = spinor_wave_functions(lapw_global.size_wfmt + ig, ispn, j1);
                
                lapw_fft(1, &v1[ispn][0]);
            }
        }
       
        for (int j2 = 0; j2 < lapw_global.nstsv; j2++)
        {
            complex16 zsum(0,0);
            
            for (int ispn = 0; ispn < lapw_global.nspinor; ispn++)
            {
                for (int ias = 0; ias < (int)lapw_global.atoms.size(); ias++)
                {
                    int offset_wfmt = lapw_global.atoms[ias]->offset_wfmt;
                    Species *sp = lapw_global.atoms[ias]->species;

                    for (int l = 0; l <= lapw_global.lmaxapw; l++)
                    {
                        int ordmax = sp->radial_index.nrf(l);
                        for (int io1 = 0; io1 < ordmax; io1++)
                            for (int io2 = 0; io2 < ordmax; io2++)
                                for (int m = -l; m <= (int)l; m++)
                                    zsum += conj(spinor_wave_functions(offset_wfmt + sp->index(l, m, io1), ispn, j1)) *
                                                 spinor_wave_functions(offset_wfmt + sp->index(l, m, io2), ispn, j2) * 
                                                 lapw_runtime.ovlprad(l, io1, io2, ias);
                    }
                }
            }

            if (use_fft == 0) 
            {
                for (int ispn = 0; ispn < lapw_global.nspinor; ispn++)
                {
                    for (int ig1 = 0; ig1 < ngk; ig1++)
                    {
                        for (int ig2 = 0; ig2 < ngk; ig2++)
                        {
                            int ig3 = idxG12(idxg[ig1], idxg[ig2]);
                            zsum += conj(spinor_wave_functions(lapw_global.size_wfmt + ig1, ispn, j1)) * spinor_wave_functions(lapw_global.size_wfmt + ig2, ispn, j2) * lapw_global.cfunig[ig3];
                        }
                    }    
                }   
            }
            if (use_fft == 1)
            {
                for (int ispn = 0; ispn < lapw_global.nspinor; ispn++)
                {
                    for (int ig = 0; ig < ngk; ig++)
                        zsum += conj(v1[ispn][ig]) * spinor_wave_functions(lapw_global.size_wfmt + ig, ispn, j2);
                }
            }
           
            if (use_fft == 2)
            {
                for (int ispn = 0; ispn < lapw_global.nspinor; ispn++)
                {
                    memset(&zfft[0], 0, lapw_global.ngrtot * sizeof(complex16));
                    for (int ig = 0; ig < ngk; ig++) 
                        zfft[idxgfft[ig]] = spinor_wave_functions(lapw_global.size_wfmt + ig, ispn, j2);
            
                    lapw_fft(1, &zfft[0]);

                    for (int ir = 0; ir < lapw_global.ngrtot; ir++)
                        zsum += conj(v1[ispn][ir]) * zfft[ir] * lapw_global.cfunir[ir] / double(lapw_global.ngrtot);
                }
            }

            zsum = (j1 == j2) ? zsum - zone : zsum;
            maxerr = std::max(maxerr, abs(zsum));
        }
    }
    std :: cout << "maximum error = " << maxerr << std::endl;
}
