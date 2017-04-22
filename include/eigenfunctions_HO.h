/* * ATUS2 - The ATUS2 package is atom interferometer Toolbox developed at ZARM
 * (CENTER OF APPLIED SPACE TECHNOLOGY AND MICROGRAVITY), Germany. This project is
 * founded by the DLR Agentur (Deutsche Luft und Raumfahrt Agentur). Grant numbers:
 * 50WM0942, 50WM1042, 50WM1342.
 * Copyright (C) 2017 Želimir Marojević, Ertan Göklü, Claus Lämmerzahl
 *
 * This file is part of ATUS2.
 *
 * ATUS2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ATUS2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ATUS2.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __eigenfunctions_HO_h__
#define __eigenfunctions_HO_h__

double HO_0( const double L, const double x )
{
  double retval = 0.7511255444e0 * pow(L, -0.1e1 / 0.2e1) * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1));
  return retval;
}

double HO_1( const double L, const double x )
{
  double retval = 0.1062251932e1 * x * pow(L, -0.3e1 / 0.2e1) * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1));
  return retval;
}

double HO_2( const double L, const double x )
{
  double retval = (-0.7511255442e0 * L * L + 0.1502251088e1 * x * x) * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.5e1 / 0.2e1);
  return retval;
}

double HO_3( const double L, const double x )
{
  double retval = (-0.3186755793e1 * x * L * L + 0.2124503862e1 * pow(x, 0.3e1)) * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.7e1 / 0.2e1);
  return retval;
}

double HO_4( const double L, const double x )
{
  double retval = exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.9e1 / 0.2e1) * (-0.9013506523e1 * x * x * L * L + 0.3004502174e1 * pow(x, 0.4e1) + 0.2253376631e1 * pow(L, 0.4e1));
  return retval;
}

double HO_5( const double L, const double x )
{
  double retval = x * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.11e2 / 0.2e1) * (-0.2124503861e2 * x * x * L * L + 0.4249007722e1 * pow(x, 0.4e1) + 0.1593377897e2 * pow(L, 0.4e1));
  return retval;
}

double HO_6( const double L, const double x )
{
  double retval = -exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.13e2 / 0.2e1) * (0.4506753258e2 * pow(x, 0.4e1) * L * L - 0.6009004346e1 * pow(x, 0.6e1) - 0.6760129891e2 * x * x * pow(L, 0.4e1) + 0.1126688316e2 * pow(L, 0.6e1));
  return retval;
}

double HO_7( const double L, const double x )
{
  double retval = -x * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.15e2 / 0.2e1) * (0.8922916210e2 * pow(x, 0.4e1) * L * L - 0.8498015440e1 * pow(x, 0.6e1) - 0.2230729053e3 * x * x * pow(L, 0.4e1) + 0.1115364527e3 * pow(L, 0.6e1));
  return retval;
}

double HO_8( const double L, const double x )
{
  double retval = exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.17e2 / 0.2e1) * (-0.1682521216e3 * pow(x, 0.6e1) * L * L + 0.1201800869e2 * pow(x, 0.8e1) + 0.6309454559e3 * pow(x, 0.4e1) * pow(L, 0.4e1) - 0.6309454561e3 * x * x * pow(L, 0.6e1) + 0.7886818203e2 * pow(L, 0.8e1));
  return retval;
}

double HO_9( const double L, const double x )
{
  double retval = x * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.19e2 / 0.2e1) * (-0.3059285557e3 * pow(x, 0.6e1) * L * L + 0.1699603088e2 * pow(x, 0.8e1) + 0.1606124917e4 * pow(x, 0.4e1) * pow(L, 0.4e1) - 0.2676874862e4 * x * x * pow(L, 0.6e1) + 0.1003828073e4 * pow(L, 0.8e1));
  return retval;
}

double HO_10( const double L, const double x )
{
  double retval = -0.7098136374e3 * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.21e2 / 0.2e1) * (0.7619047622e0 * pow(x, 0.8e1) * L * L - 0.3386243390e-1 * pow(x, 0.10e2) - 0.5333333335e1 * pow(x, 0.6e1) * pow(L, 0.4e1) + 0.1333333334e2 * pow(x, 0.4e1) * pow(L, 0.6e1) - 0.1000000000e2 * x * x * pow(L, 0.8e1) + 0.1000000000e1 * pow(L, 0.10e2));
  return retval;
}

double HO_11( const double L, const double x )
{
  double retval = -0.1104210881e5 * x * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.23e2 / 0.2e1) * (0.8465608462e-1 * pow(x, 0.8e1) * L * L - 0.3078403079e-2 * pow(x, 0.10e2) - 0.7619047615e0 * pow(x, 0.6e1) * pow(L, 0.4e1) + 0.2666666665e1 * pow(x, 0.4e1) * pow(L, 0.6e1) - 0.3333333331e1 * x * x * pow(L, 0.8e1) + 0.1000000000e1 * pow(L, 0.10e2));
  return retval;
}

double HO_12( const double L, const double x )
{
  double retval = 0.7807950016e4 * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.25e2 / 0.2e1) * (-0.2031746031e0 * pow(x, 0.10e2) * L * L + 0.6156806159e-2 * pow(x, 0.12e2) + 0.2285714286e1 * pow(x, 0.8e1) * pow(L, 0.4e1) - 0.1066666666e2 * pow(x, 0.6e1) * pow(L, 0.6e1) + 0.2000000000e2 * pow(x, 0.4e1) * pow(L, 0.8e1) - 0.1199999999e2 * x * x * pow(L, 0.10e2) + 0.1000000000e1 * pow(L, 0.12e2));
  return retval;
}

double HO_13( const double L, const double x )
{
  double retval = 0.1435474143e6 * x * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.27e2 / 0.2e1) * (-0.1847041848e-1 * pow(x, 0.10e2) * L * L + 0.4736004742e-3 * pow(x, 0.12e2) + 0.2539682542e0 * pow(x, 0.8e1) * pow(L, 0.4e1) - 0.1523809525e1 * pow(x, 0.6e1) * pow(L, 0.6e1) + 0.4000000002e1 * pow(x, 0.4e1) * pow(L, 0.8e1) - 0.4000000003e1 * x * x * pow(L, 0.10e2) + 0.1000000000e1 * pow(L, 0.12e2));
  return retval;
}

double HO_14( const double L, const double x )
{
  double retval = -0.1015033500e6 * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.29e2 / 0.2e1) * (0.4309764314e-1 * pow(x, 0.12e2) * L * L - 0.9472009486e-3 * pow(x, 0.14e2) - 0.7111111118e0 * pow(x, 0.10e2) * pow(L, 0.4e1) + 0.5333333339e1 * pow(x, 0.8e1) * pow(L, 0.6e1) - 0.1866666669e2 * pow(x, 0.6e1) * pow(L, 0.8e1) + 0.2800000003e2 * pow(x, 0.4e1) * pow(L, 0.10e2) - 0.1400000002e2 * x * x * pow(L, 0.12e2) + 0.1000000000e1 * pow(L, 0.14e2));
  return retval;
}

double HO_15( const double L, const double x )
{
  double retval = -0.2153211215e7 * x * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.31e2 / 0.2e1) * (0.3315203315e-2 * pow(x, 0.12e2) * L * L - 0.6314672984e-4 * pow(x, 0.14e2) - 0.6464646465e-1 * pow(x, 0.10e2) * pow(L, 0.4e1) + 0.5925925924e0 * pow(x, 0.8e1) * pow(L, 0.6e1) - 0.2666666666e1 * pow(x, 0.6e1) * pow(L, 0.8e1) + 0.5599999998e1 * pow(x, 0.4e1) * pow(L, 0.10e2) - 0.4666666665e1 * x * x * pow(L, 0.12e2) + 0.9999999999e0 * pow(L, 0.14e2));
  return retval;
}

double HO_16( const double L, const double x )
{
  double retval = 0.1522550251e7 * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.33e2 / 0.2e1) * (-0.7577607578e-2 * pow(x, 0.14e2) * L * L + 0.1262934597e-3 * pow(x, 0.16e2) + 0.1723905724e0 * pow(x, 0.12e2) * pow(L, 0.4e1) - 0.1896296296e1 * pow(x, 0.10e2) * pow(L, 0.6e1) + 0.1066666667e2 * pow(x, 0.8e1) * pow(L, 0.8e1) - 0.2986666667e2 * pow(x, 0.6e1) * pow(L, 0.10e2) + 0.3733333332e2 * pow(x, 0.4e1) * pow(L, 0.12e2) - 0.1600000000e2 * x * x * pow(L, 0.14e2) + 0.1000000000e1 * pow(L, 0.16e2));
  return retval;
}

double HO_17( const double L, const double x )
{
  double retval = 0.3660459064e8 * x * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.35e2 / 0.2e1) * (-0.5051738386e-3 * pow(x, 0.14e2) * L * L + 0.7429027043e-5 * pow(x, 0.16e2) + 0.1326081326e-1 * pow(x, 0.12e2) * pow(L, 0.4e1) - 0.1723905724e0 * pow(x, 0.10e2) * pow(L, 0.6e1) + 0.1185185185e1 * pow(x, 0.8e1) * pow(L, 0.8e1) - 0.4266666667e1 * pow(x, 0.6e1) * pow(L, 0.10e2) + 0.7466666669e1 * pow(x, 0.4e1) * pow(L, 0.12e2) - 0.5333333331e1 * x * x * pow(L, 0.14e2) + 0.1000000000e1 * pow(L, 0.16e2));
  return retval;
}

double HO_18( const double L, const double x )
{
  double retval = -0.2588335426e8 * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.37e2 / 0.2e1) * (-0.1485805409e-4 * pow(x, 0.18e2) + 0.1136641137e-2 * pow(x, 0.16e2) * L * L - 0.3409923410e-1 * pow(x, 0.14e2) * pow(L, 0.4e1) + 0.5171717173e0 * pow(x, 0.12e2) * pow(L, 0.6e1) - 0.4266666665e1 * pow(x, 0.10e2) * pow(L, 0.8e1) + 0.1920000000e2 * pow(x, 0.8e1) * pow(L, 0.10e2) - 0.4480000001e2 * pow(x, 0.6e1) * pow(L, 0.12e2) + 0.4800000002e2 * pow(x, 0.4e1) * pow(L, 0.14e2) - 0.1799999999e2 * x * x * pow(L, 0.16e2) + 0.1000000000e1 * pow(L, 0.18e2));
  return retval;
}

double HO_19( const double L, const double x )
{
  double retval = -0.6954872214e9 * x * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.39e2 / 0.2e1) * (-0.7820028474e-6 * pow(x, 0.18e2) + 0.6686124342e-4 * pow(x, 0.16e2) * L * L - 0.2273282276e-2 * pow(x, 0.14e2) * pow(L, 0.4e1) + 0.3978243983e-1 * pow(x, 0.12e2) * pow(L, 0.6e1) - 0.3878787882e0 * pow(x, 0.10e2) * pow(L, 0.8e1) + 0.2133333335e1 * pow(x, 0.8e1) * pow(L, 0.10e2) - 0.6400000005e1 * pow(x, 0.6e1) * pow(L, 0.12e2) + 0.9600000011e1 * pow(x, 0.4e1) * pow(L, 0.14e2) - 0.6000000006e1 * x * x * pow(L, 0.16e2) + 0.1000000000e1 * pow(L, 0.18e2));
  return retval;
}

double HO_20( const double L, const double x )
{
  double retval = 0.4917837304e9 * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.41e2 / 0.2e1) * (0.1564005695e-5 * pow(x, 0.20e2) - 0.1485805410e-3 * pow(x, 0.18e2) * L * L + 0.5683205691e-2 * pow(x, 0.16e2) * pow(L, 0.4e1) - 0.1136641138e0 * pow(x, 0.14e2) * pow(L, 0.6e1) + 0.1292929294e1 * pow(x, 0.12e2) * pow(L, 0.8e1) - 0.8533333339e1 * pow(x, 0.10e2) * pow(L, 0.10e2) + 0.3200000004e2 * pow(x, 0.8e1) * pow(L, 0.12e2) - 0.6400000005e2 * pow(x, 0.6e1) * pow(L, 0.14e2) + 0.6000000004e2 * pow(x, 0.4e1) * pow(L, 0.16e2) - 0.2000000002e2 * x * x * pow(L, 0.18e2) + 0.1000000000e1 * pow(L, 0.20e2));
  return retval;
}

double HO_21( const double L, const double x )
{
  double retval = 0.1460523166e11 * x * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.43e2 / 0.2e1) * (0.7447646161e-7 * pow(x, 0.20e2) - 0.7820028464e-5 * pow(x, 0.18e2) * L * L + 0.3343062168e-3 * pow(x, 0.16e2) * pow(L, 0.4e1) - 0.7577607578e-2 * pow(x, 0.14e2) * pow(L, 0.6e1) + 0.9945609942e-1 * pow(x, 0.12e2) * pow(L, 0.8e1) - 0.7757575754e0 * pow(x, 0.10e2) * pow(L, 0.10e2) + 0.3555555554e1 * pow(x, 0.8e1) * pow(L, 0.12e2) - 0.9142857143e1 * pow(x, 0.6e1) * pow(L, 0.14e2) + 0.1199999999e2 * pow(x, 0.4e1) * pow(L, 0.16e2) - 0.6666666666e1 * x * x * pow(L, 0.18e2) + 0.1000000000e1 * pow(L, 0.20e2));
  return retval;
}

double HO_22( const double L, const double x )
{
  double retval = -0.1032745834e11 * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.45e2 / 0.2e1) * (0.1720406262e-4 * pow(x, 0.20e2) * L * L - 0.8171929747e-3 * pow(x, 0.18e2) * pow(L, 0.4e1) + 0.2083842085e-1 * pow(x, 0.16e2) * pow(L, 0.6e1) - 0.3125763127e0 * pow(x, 0.14e2) * pow(L, 0.8e1) + 0.2844444444e1 * pow(x, 0.12e2) * pow(L, 0.10e2) - 0.1564444445e2 * pow(x, 0.10e2) * pow(L, 0.12e2) + 0.5028571428e2 * pow(x, 0.8e1) * pow(L, 0.14e2) - 0.8800000003e2 * pow(x, 0.6e1) * pow(L, 0.16e2) + 0.7333333332e2 * pow(x, 0.4e1) * pow(L, 0.18e2) - 0.2200000000e2 * x * x * pow(L, 0.20e2) - 0.1489529233e-6 * pow(x, 0.22e2) + 0.1000000000e1 * pow(L, 0.22e2));
  return retval;
}

double HO_23( const double L, const double x )
{
  double retval = -0.3359203278e12 * x * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.47e2 / 0.2e1) * (0.8192410773e-6 * pow(x, 0.20e2) * L * L - 0.4301015656e-4 * pow(x, 0.18e2) * pow(L, 0.4e1) + 0.1225789462e-2 * pow(x, 0.16e2) * pow(L, 0.6e1) - 0.2083842086e-1 * pow(x, 0.14e2) * pow(L, 0.8e1) + 0.2188034189e0 * pow(x, 0.12e2) * pow(L, 0.10e2) - 0.1422222223e1 * pow(x, 0.10e2) * pow(L, 0.12e2) + 0.5587301591e1 * pow(x, 0.8e1) * pow(L, 0.14e2) - 0.1257142858e2 * pow(x, 0.6e1) * pow(L, 0.16e2) + 0.1466666667e2 * pow(x, 0.4e1) * pow(L, 0.18e2) - 0.7333333333e1 * x * x * pow(L, 0.20e2) - 0.6476214061e-8 * pow(x, 0.22e2) + 0.1000000000e1 * pow(L, 0.22e2));
  return retval;
}

double HO_24( const double L, const double x )
{
  double retval = 0.2375315417e12 * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.49e2 / 0.2e1) * (0.1295242812e-7 * pow(x, 0.24e2) - 0.1787435077e-5 * pow(x, 0.22e2) * L * L + 0.1032243757e-3 * pow(x, 0.20e2) * pow(L, 0.4e1) - 0.3268771898e-2 * pow(x, 0.18e2) * pow(L, 0.6e1) + 0.6251526258e-1 * pow(x, 0.16e2) * pow(L, 0.8e1) - 0.7501831505e0 * pow(x, 0.14e2) * pow(L, 0.10e2) + 0.5688888892e1 * pow(x, 0.12e2) * pow(L, 0.12e2) - 0.2681904763e2 * pow(x, 0.10e2) * pow(L, 0.14e2) + 0.7542857147e2 * pow(x, 0.8e1) * pow(L, 0.16e2) - 0.1173333334e3 * pow(x, 0.6e1) * pow(L, 0.18e2) + 0.8800000000e2 * pow(x, 0.4e1) * pow(L, 0.20e2) - 0.2400000000e2 * x * x * pow(L, 0.22e2) + 0.1000000000e1 * pow(L, 0.24e2));
  return retval;
}

double HO_25( const double L, const double x )
{
  double retval = 0.8398008189e13 * x * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.51e2 / 0.2e1) * (0.5180971247e-9 * pow(x, 0.24e2) - 0.7771456856e-7 * pow(x, 0.22e2) * L * L + 0.4915446460e-5 * pow(x, 0.20e2) * pow(L, 0.4e1) - 0.1720406262e-3 * pow(x, 0.18e2) * pow(L, 0.6e1) + 0.3677368384e-2 * pow(x, 0.16e2) * pow(L, 0.8e1) - 0.5001221002e-1 * pow(x, 0.14e2) * pow(L, 0.10e2) + 0.4376068378e0 * pow(x, 0.12e2) * pow(L, 0.12e2) - 0.2438095240e1 * pow(x, 0.10e2) * pow(L, 0.14e2) + 0.8380952383e1 * pow(x, 0.8e1) * pow(L, 0.16e2) - 0.1676190477e2 * pow(x, 0.6e1) * pow(L, 0.18e2) + 0.1760000000e2 * pow(x, 0.4e1) * pow(L, 0.20e2) - 0.7999999998e1 * x * x * pow(L, 0.22e2) + 0.9999999998e0 * pow(L, 0.24e2));
  return retval;
}

double HO_26( const double L, const double x )
{
  double retval = -0.5938288536e13 * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.53e2 / 0.2e1) * (0.1000000000e1 * pow(L, 0.26e2) - 0.1036194250e-8 * pow(x, 0.26e2) + 0.1683815653e-6 * pow(x, 0.24e2) * L * L - 0.1161832800e-4 * pow(x, 0.22e2) * pow(L, 0.4e1) + 0.4473056283e-3 * pow(x, 0.20e2) * pow(L, 0.6e1) - 0.1062350867e-1 * pow(x, 0.18e2) * pow(L, 0.8e1) + 0.1625396826e0 * pow(x, 0.16e2) * pow(L, 0.10e2) - 0.1625396826e1 * pow(x, 0.14e2) * pow(L, 0.12e2) + 0.1056507938e2 * pow(x, 0.12e2) * pow(L, 0.14e2) - 0.4358095244e2 * pow(x, 0.10e2) * pow(L, 0.16e2) + 0.1089523810e3 * pow(x, 0.8e1) * pow(L, 0.18e2) - 0.1525333334e3 * pow(x, 0.6e1) * pow(L, 0.20e2) + 0.1040000000e3 * pow(x, 0.4e1) * pow(L, 0.22e2) - 0.2600000000e2 * x * x * pow(L, 0.24e2));
  return retval;
}

double HO_27( const double L, const double x )
{
  double retval = -0.2267462210e15 * x * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.55e2 / 0.2e1) * (0.1000000000e1 * pow(L, 0.26e2) - 0.3837756481e-10 * pow(x, 0.26e2) + 0.6735262609e-8 * pow(x, 0.24e2) * L * L - 0.5051446957e-6 * pow(x, 0.22e2) * pow(L, 0.4e1) + 0.2130026801e-4 * pow(x, 0.20e2) * pow(L, 0.6e1) - 0.5591320351e-3 * pow(x, 0.18e2) * pow(L, 0.8e1) + 0.9561157798e-2 * pow(x, 0.16e2) * pow(L, 0.10e2) - 0.1083597883e0 * pow(x, 0.14e2) * pow(L, 0.12e2) + 0.8126984132e0 * pow(x, 0.12e2) * pow(L, 0.14e2) - 0.3961904767e1 * pow(x, 0.10e2) * pow(L, 0.16e2) + 0.1210582012e2 * pow(x, 0.8e1) * pow(L, 0.18e2) - 0.2179047619e2 * pow(x, 0.6e1) * pow(L, 0.20e2) + 0.2080000001e2 * pow(x, 0.4e1) * pow(L, 0.22e2) - 0.8666666661e1 * x * x * pow(L, 0.24e2));
  return retval;
}

double HO_28( const double L, const double x )
{
  double retval = 0.1603337904e15 * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.57e2 / 0.2e1) * (-0.2799999999e2 * x * x * pow(L, 0.26e2) - 0.1450671947e-7 * pow(x, 0.26e2) * L * L + 0.1178670957e-5 * pow(x, 0.24e2) * pow(L, 0.4e1) - 0.5421886404e-4 * pow(x, 0.22e2) * pow(L, 0.6e1) + 0.1565569699e-2 * pow(x, 0.20e2) * pow(L, 0.8e1) - 0.2974582426e-1 * pow(x, 0.18e2) * pow(L, 0.10e2) + 0.3792592593e0 * pow(x, 0.16e2) * pow(L, 0.12e2) - 0.3250793651e1 * pow(x, 0.14e2) * pow(L, 0.14e2) + 0.1848888891e2 * pow(x, 0.12e2) * pow(L, 0.16e2) - 0.6779259271e2 * pow(x, 0.10e2) * pow(L, 0.18e2) + 0.1525333335e3 * pow(x, 0.8e1) * pow(L, 0.20e2) - 0.1941333334e3 * pow(x, 0.6e1) * pow(L, 0.22e2) + 0.1213333334e3 * pow(x, 0.4e1) * pow(L, 0.24e2) + 0.7675512965e-10 * pow(x, 0.28e2) + 0.1000000000e1 * pow(L, 0.28e2));
  return retval;
}

double HO_29( const double L, const double x )
{
  double retval = 0.6575640401e16 * x * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.59e2 / 0.2e1) * (-0.9333333342e1 * x * x * pow(L, 0.26e2) - 0.5372859066e-9 * pow(x, 0.26e2) * L * L + 0.4714683830e-7 * pow(x, 0.24e2) * pow(L, 0.4e1) - 0.2357341916e-5 * pow(x, 0.22e2) * pow(L, 0.6e1) + 0.7455093808e-4 * pow(x, 0.20e2) * pow(L, 0.8e1) - 0.1565569698e-2 * pow(x, 0.18e2) * pow(L, 0.10e2) + 0.2230936821e-1 * pow(x, 0.16e2) * pow(L, 0.12e2) - 0.2167195768e0 * pow(x, 0.14e2) * pow(L, 0.14e2) + 0.1422222224e1 * pow(x, 0.12e2) * pow(L, 0.16e2) - 0.6162962975e1 * pow(x, 0.10e2) * pow(L, 0.18e2) + 0.1694814818e2 * pow(x, 0.8e1) * pow(L, 0.20e2) - 0.2773333338e2 * pow(x, 0.6e1) * pow(L, 0.22e2) + 0.2426666668e2 * pow(x, 0.4e1) * pow(L, 0.24e2) + 0.2646728610e-11 * pow(x, 0.28e2) + 0.1000000000e1 * pow(L, 0.28e2));
  return retval;
}

double HO_30( const double L, const double x )
{
  double retval = -0.4649679917e16 * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.61e2 / 0.2e1) * (0.1400000001e3 * pow(x, 0.4e1) * pow(L, 0.26e2) + 0.1151326943e-8 * pow(x, 0.28e2) * L * L - 0.1088003961e-6 * pow(x, 0.26e2) * pow(L, 0.4e1) + 0.5893354788e-5 * pow(x, 0.24e2) * pow(L, 0.6e1) - 0.2033207402e-3 * pow(x, 0.22e2) * pow(L, 0.8e1) + 0.4696709096e-2 * pow(x, 0.20e2) * pow(L, 0.10e2) - 0.7436456067e-1 * pow(x, 0.18e2) * pow(L, 0.12e2) + 0.8126984133e0 * pow(x, 0.16e2) * pow(L, 0.14e2) - 0.6095238099e1 * pow(x, 0.14e2) * pow(L, 0.16e2) + 0.3081481486e2 * pow(x, 0.12e2) * pow(L, 0.18e2) - 0.1016888891e3 * pow(x, 0.10e2) * pow(L, 0.20e2) + 0.2080000004e3 * pow(x, 0.8e1) * pow(L, 0.22e2) - 0.2426666670e3 * pow(x, 0.6e1) * pow(L, 0.24e2) - 0.3000000002e2 * x * x * pow(L, 0.28e2) - 0.5293457218e-11 * pow(x, 0.30e2) + 0.1000000000e1 * pow(L, 0.30e2));
  return retval;
}

double HO_31( const double L, const double x )
{
  double retval = -0.2038448525e18 * x * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.63e2 / 0.2e1) * (0.1000000000e1 * pow(L, 0.30e2) + 0.2800000001e2 * pow(x, 0.4e1) * pow(L, 0.26e2) + 0.3970092904e-10 * pow(x, 0.28e2) * L * L - 0.4029644297e-8 * pow(x, 0.26e2) * pow(L, 0.4e1) + 0.2357341913e-6 * pow(x, 0.24e2) * pow(L, 0.6e1) - 0.8840032176e-5 * pow(x, 0.22e2) * pow(L, 0.8e1) + 0.2236528139e-3 * pow(x, 0.20e2) * pow(L, 0.10e2) - 0.3913924243e-2 * pow(x, 0.18e2) * pow(L, 0.12e2) + 0.4780578897e-1 * pow(x, 0.16e2) * pow(L, 0.14e2) - 0.4063492062e0 * pow(x, 0.14e2) * pow(L, 0.16e2) + 0.2370370370e1 * pow(x, 0.12e2) * pow(L, 0.18e2) - 0.9244444453e1 * pow(x, 0.10e2) * pow(L, 0.20e2) + 0.2311111114e2 * pow(x, 0.8e1) * pow(L, 0.22e2) - 0.3466666670e2 * pow(x, 0.6e1) * pow(L, 0.24e2) - 0.9999999995e1 * x * x * pow(L, 0.28e2) - 0.1707566843e-12 * pow(x, 0.30e2));
  return retval;
}

double HO_32( const double L, const double x )
{
  double retval = 0.1441400775e18 * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.65e2 / 0.2e1) * (0.3415133684e-12 * pow(x, 0.32e2) - 0.3199999998e2 * x * x * pow(L, 0.30e2) - 0.2986666669e3 * pow(x, 0.6e1) * pow(L, 0.26e2) - 0.8469531523e-10 * pow(x, 0.30e2) * L * L + 0.9210615535e-8 * pow(x, 0.28e2) * pow(L, 0.4e1) - 0.5802687784e-6 * pow(x, 0.26e2) * pow(L, 0.6e1) + 0.2357341913e-4 * pow(x, 0.24e2) * pow(L, 0.8e1) - 0.6506263677e-3 * pow(x, 0.22e2) * pow(L, 0.10e2) + 0.1252455757e-1 * pow(x, 0.20e2) * pow(L, 0.12e2) - 0.1699761385e0 * pow(x, 0.18e2) * pow(L, 0.14e2) + 0.1625396824e1 * pow(x, 0.16e2) * pow(L, 0.16e2) - 0.1083597883e2 * pow(x, 0.14e2) * pow(L, 0.18e2) + 0.4930370371e2 * pow(x, 0.12e2) * pow(L, 0.20e2) - 0.1479111113e3 * pow(x, 0.10e2) * pow(L, 0.22e2) + 0.2773333336e3 * pow(x, 0.8e1) * pow(L, 0.24e2) + 0.1600000000e3 * pow(x, 0.4e1) * pow(L, 0.28e2) + 0.1000000000e1 * pow(L, 0.32e2));
  return retval;
}

double HO_33( const double L, const double x )
{
  double retval = 0.6726880125e19 * x * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.67e2 / 0.2e1) * (0.1034888996e-13 * pow(x, 0.32e2) - 0.1066666668e2 * x * x * pow(L, 0.30e2) - 0.4266666675e2 * pow(x, 0.6e1) * pow(L, 0.26e2) - 0.2732106946e-11 * pow(x, 0.30e2) * L * L + 0.3176074326e-9 * pow(x, 0.28e2) * pow(L, 0.4e1) - 0.2149143626e-7 * pow(x, 0.26e2) * pow(L, 0.6e1) + 0.9429367662e-6 * pow(x, 0.24e2) * pow(L, 0.8e1) - 0.2828810297e-4 * pow(x, 0.22e2) * pow(L, 0.10e2) + 0.5964075042e-3 * pow(x, 0.20e2) * pow(L, 0.12e2) - 0.8946112562e-2 * pow(x, 0.18e2) * pow(L, 0.14e2) + 0.9561157799e-1 * pow(x, 0.16e2) * pow(L, 0.16e2) - 0.7223985892e0 * pow(x, 0.14e2) * pow(L, 0.18e2) + 0.3792592595e1 * pow(x, 0.12e2) * pow(L, 0.20e2) - 0.1344646467e2 * pow(x, 0.10e2) * pow(L, 0.22e2) + 0.3081481489e2 * pow(x, 0.8e1) * pow(L, 0.24e2) + 0.3200000005e2 * pow(x, 0.4e1) * pow(L, 0.28e2) + 0.1000000000e1 * pow(L, 0.32e2));
  return retval;
}

double HO_34( const double L, const double x )
{
  double retval = -0.4756622551e19 * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.69e2 / 0.2e1) * (0.1813333337e3 * pow(x, 0.4e1) * pow(L, 0.30e2) + 0.3626666677e3 * pow(x, 0.8e1) * pow(L, 0.26e2) + 0.5805727260e-11 * pow(x, 0.32e2) * L * L - 0.7199101807e-9 * pow(x, 0.30e2) * pow(L, 0.4e1) + 0.5219348808e-7 * pow(x, 0.28e2) * pow(L, 0.6e1) - 0.2466142311e-5 * pow(x, 0.26e2) * pow(L, 0.8e1) + 0.8014962510e-4 * pow(x, 0.24e2) * pow(L, 0.10e2) - 0.1843441377e-2 * pow(x, 0.22e2) * pow(L, 0.12e2) + 0.3041678272e-1 * pow(x, 0.20e2) * pow(L, 0.14e2) - 0.3611992946e0 * pow(x, 0.18e2) * pow(L, 0.16e2) + 0.3070194005e1 * pow(x, 0.16e2) * pow(L, 0.18e2) - 0.1842116403e2 * pow(x, 0.14e2) * pow(L, 0.20e2) + 0.7619663310e2 * pow(x, 0.12e2) * pow(L, 0.22e2) - 0.2095407412e3 * pow(x, 0.10e2) * pow(L, 0.24e2) - 0.3626666675e3 * pow(x, 0.6e1) * pow(L, 0.28e2) - 0.3400000006e2 * x * x * pow(L, 0.32e2) - 0.2069777992e-13 * pow(x, 0.34e2) + 0.1000000000e1 * pow(L, 0.34e2));
  return retval;
}

double HO_35( const double L, const double x )
{
  double retval = -0.2354408047e21 * x * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.71e2 / 0.2e1) * (0.3626666666e2 * pow(x, 0.4e1) * pow(L, 0.30e2) + 0.4029629631e2 * pow(x, 0.8e1) * pow(L, 0.26e2) + 0.1759311287e-12 * pow(x, 0.32e2) * L * L - 0.2322290901e-10 * pow(x, 0.30e2) * pow(L, 0.4e1) + 0.1799775447e-8 * pow(x, 0.28e2) * pow(L, 0.6e1) - 0.9133860389e-7 * pow(x, 0.26e2) * pow(L, 0.8e1) + 0.3205984998e-5 * pow(x, 0.24e2) * pow(L, 0.10e2) - 0.8014962493e-4 * pow(x, 0.22e2) * pow(L, 0.12e2) + 0.1448418221e-2 * pow(x, 0.20e2) * pow(L, 0.14e2) - 0.1901048915e-1 * pow(x, 0.18e2) * pow(L, 0.16e2) + 0.1805996470e0 * pow(x, 0.16e2) * pow(L, 0.18e2) - 0.1228077599e1 * pow(x, 0.14e2) * pow(L, 0.20e2) + 0.5861279453e1 * pow(x, 0.12e2) * pow(L, 0.22e2) - 0.1904915824e2 * pow(x, 0.10e2) * pow(L, 0.24e2) - 0.5180952382e2 * pow(x, 0.6e1) * pow(L, 0.28e2) - 0.1133333333e2 * x * x * pow(L, 0.32e2) - 0.5913651395e-15 * pow(x, 0.34e2) + 0.1000000000e1 * pow(L, 0.34e2));
  return retval;
}

double HO_36( const double L, const double x )
{
  double retval = 0.1664817895e21 * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.73e2 / 0.2e1) * (0.6843776094e-1 * pow(x, 0.20e2) * pow(L, 0.16e2) - 0.7223985882e0 * pow(x, 0.18e2) * pow(L, 0.18e2) + 0.5526349198e1 * pow(x, 0.16e2) * pow(L, 0.20e2) - 0.3014372290e2 * pow(x, 0.14e2) * pow(L, 0.22e2) + 0.1142949494e3 * pow(x, 0.12e2) * pow(L, 0.24e2) + 0.4662857145e3 * pow(x, 0.8e1) * pow(L, 0.28e2) + 0.2040000000e3 * pow(x, 0.4e1) * pow(L, 0.32e2) - 0.3600000000e2 * x * x * pow(L, 0.34e2) + 0.1182730279e-14 * pow(x, 0.36e2) + 0.1000000000e1 * pow(L, 0.36e2) - 0.4352000002e3 * pow(x, 0.6e1) * pow(L, 0.30e2) - 0.2901333334e3 * pow(x, 0.10e2) * pow(L, 0.26e2) - 0.3725600374e-12 * pow(x, 0.34e2) * L * L + 0.5225154529e-10 * pow(x, 0.32e2) * pow(L, 0.4e1) - 0.4319461074e-8 * pow(x, 0.30e2) * pow(L, 0.6e1) + 0.2348706959e-6 * pow(x, 0.28e2) * pow(L, 0.8e1) - 0.8878112306e-5 * pow(x, 0.26e2) * pow(L, 0.10e2) + 0.2404488749e-3 * pow(x, 0.24e2) * pow(L, 0.12e2) - 0.4740277816e-2 * pow(x, 0.22e2) * pow(L, 0.14e2));
  return retval;
}

double HO_37( const double L, const double x )
{
  double retval = 0.8711309762e22 * x * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.75e2 / 0.2e1) * (0.3258940998e-2 * pow(x, 0.20e2) * pow(L, 0.16e2) - 0.3802097831e-1 * pow(x, 0.18e2) * pow(L, 0.18e2) + 0.3250793647e0 * pow(x, 0.16e2) * pow(L, 0.20e2) - 0.2009581526e1 * pow(x, 0.14e2) * pow(L, 0.22e2) + 0.8791919179e1 * pow(x, 0.12e2) * pow(L, 0.24e2) + 0.5180952381e2 * pow(x, 0.8e1) * pow(L, 0.28e2) + 0.4080000001e2 * pow(x, 0.4e1) * pow(L, 0.32e2) - 0.1200000000e2 * x * x * pow(L, 0.34e2) + 0.3196568320e-16 * pow(x, 0.36e2) + 0.9999999996e0 * pow(L, 0.36e2) - 0.6217142859e2 * pow(x, 0.6e1) * pow(L, 0.30e2) - 0.2637575756e2 * pow(x, 0.10e2) * pow(L, 0.26e2) - 0.1064457249e-13 * pow(x, 0.34e2) * L * L + 0.1583380159e-11 * pow(x, 0.32e2) * pow(L, 0.4e1) - 0.1393374540e-9 * pow(x, 0.30e2) * pow(L, 0.6e1) + 0.8098989511e-8 * pow(x, 0.28e2) * pow(L, 0.8e1) - 0.3288189741e-6 * pow(x, 0.26e2) * pow(L, 0.10e2) + 0.9617954995e-5 * pow(x, 0.24e2) * pow(L, 0.12e2) - 0.2060990355e-3 * pow(x, 0.22e2) * pow(L, 0.14e2));
  return retval;
}

double HO_38( const double L, const double x )
{
  double retval = -0.6159826202e22 * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.77e2 / 0.2e1) * (0.1000000000e1 * pow(L, 0.38e2) - 0.1125815982e-1 * pow(x, 0.22e2) * pow(L, 0.16e2) + 0.1444797177e0 * pow(x, 0.20e2) * pow(L, 0.18e2) - 0.1372557318e1 * pow(x, 0.18e2) * pow(L, 0.20e2) + 0.9545512256e1 * pow(x, 0.16e2) * pow(L, 0.22e2) - 0.4772756127e2 * pow(x, 0.14e2) * pow(L, 0.24e2) - 0.3937523809e3 * pow(x, 0.10e2) * pow(L, 0.28e2) - 0.5168000003e3 * pow(x, 0.6e1) * pow(L, 0.32e2) + 0.2279999999e3 * pow(x, 0.4e1) * pow(L, 0.34e2) - 0.3800000001e2 * x * x * pow(L, 0.36e2) + 0.5906285718e3 * pow(x, 0.8e1) * pow(L, 0.30e2) + 0.1670464645e3 * pow(x, 0.12e2) * pow(L, 0.26e2) + 0.2247187528e-13 * pow(x, 0.36e2) * L * L - 0.3539320359e-11 * pow(x, 0.34e2) * pow(L, 0.4e1) + 0.3309264534e-9 * pow(x, 0.32e2) * pow(L, 0.6e1) - 0.2051744012e-7 * pow(x, 0.30e2) * pow(L, 0.8e1) + 0.8925086444e-6 * pow(x, 0.28e2) * pow(L, 0.10e2) - 0.2811402231e-4 * pow(x, 0.26e2) * pow(L, 0.12e2) + 0.6526469464e-3 * pow(x, 0.24e2) * pow(L, 0.14e2) - 0.6393136644e-16 * pow(x, 0.38e2));
  return retval;
}

double HO_39( const double L, const double x )
{
  double retval = -0.3397410806e24 * x * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.79e2 / 0.2e1) * (0.1000000000e1 * pow(L, 0.38e2) - 0.4894852092e-3 * pow(x, 0.22e2) * pow(L, 0.16e2) + 0.6879986554e-2 * pow(x, 0.20e2) * pow(L, 0.18e2) - 0.7223985883e-1 * pow(x, 0.18e2) * pow(L, 0.20e2) + 0.5615007204e0 * pow(x, 0.16e2) * pow(L, 0.22e2) - 0.3181837416e1 * pow(x, 0.14e2) * pow(L, 0.24e2) - 0.3579567095e2 * pow(x, 0.10e2) * pow(L, 0.28e2) - 0.7382857145e2 * pow(x, 0.6e1) * pow(L, 0.32e2) + 0.4559999999e2 * pow(x, 0.4e1) * pow(L, 0.34e2) - 0.1266666665e2 * x * x * pow(L, 0.36e2) + 0.6562539679e2 * pow(x, 0.8e1) * pow(L, 0.30e2) + 0.1284972803e2 * pow(x, 0.12e2) * pow(L, 0.26e2) + 0.6073479806e-15 * pow(x, 0.36e2) * L * L - 0.1011234387e-12 * pow(x, 0.34e2) * pow(L, 0.4e1) + 0.1002807434e-10 * pow(x, 0.32e2) * pow(L, 0.6e1) - 0.6618529067e-9 * pow(x, 0.30e2) * pow(L, 0.8e1) + 0.3077616013e-7 * pow(x, 0.28e2) * pow(L, 0.10e2) - 0.1041260085e-5 * pow(x, 0.26e2) * pow(L, 0.12e2) + 0.2610587784e-4 * pow(x, 0.24e2) * pow(L, 0.14e2) - 0.1639265805e-17 * pow(x, 0.38e2));
  return retval;
}

double HO_40( const double L, const double x )
{
  double retval = 0.2402332219e24 * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.81e2 / 0.2e1) * (0.1000000000e1 * pow(L, 0.40e2) - 0.3999999996e2 * x * x * pow(L, 0.38e2) + 0.1631617365e-2 * pow(x, 0.24e2) * pow(L, 0.16e2) - 0.2501813291e-1 * pow(x, 0.22e2) * pow(L, 0.18e2) + 0.2889594353e0 * pow(x, 0.20e2) * pow(L, 0.20e2) - 0.2495558758e1 * pow(x, 0.18e2) * pow(L, 0.22e2) + 0.1590918708e2 * pow(x, 0.16e2) * pow(L, 0.24e2) + 0.2386378063e3 * pow(x, 0.12e2) * pow(L, 0.28e2) + 0.7382857141e3 * pow(x, 0.8e1) * pow(L, 0.32e2) - 0.6079999999e3 * pow(x, 0.6e1) * pow(L, 0.34e2) + 0.2533333332e3 * pow(x, 0.4e1) * pow(L, 0.36e2) - 0.5250031740e3 * pow(x, 0.10e2) * pow(L, 0.30e2) - 0.7342701730e2 * pow(x, 0.14e2) * pow(L, 0.26e2) - 0.1278627327e-14 * pow(x, 0.38e2) * L * L + 0.2247187527e-12 * pow(x, 0.36e2) * pow(L, 0.4e1) - 0.2359546903e-10 * pow(x, 0.34e2) * pow(L, 0.6e1) + 0.1654632266e-8 * pow(x, 0.32e2) * pow(L, 0.8e1) - 0.8206976035e-7 * pow(x, 0.30e2) * pow(L, 0.10e2) + 0.2975028814e-5 * pow(x, 0.28e2) * pow(L, 0.12e2) - 0.8032577795e-4 * pow(x, 0.26e2) * pow(L, 0.14e2) + 0.3278531610e-17 * pow(x, 0.40e2));
  return retval;
}

double HO_41( const double L, const double x )
{
  double retval = 0.1392938428e26 * x * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.83e2 / 0.2e1) * (0.9999999999e0 * pow(L, 0.40e2) - 0.1333333334e2 * x * x * pow(L, 0.38e2) + 0.6526469468e-4 * pow(x, 0.24e2) * pow(L, 0.16e2) - 0.1087744911e-2 * pow(x, 0.22e2) * pow(L, 0.18e2) + 0.1375997312e-1 * pow(x, 0.20e2) * pow(L, 0.20e2) - 0.1313451980e0 * pow(x, 0.18e2) * pow(L, 0.22e2) + 0.9358345349e0 * pow(x, 0.16e2) * pow(L, 0.24e2) + 0.1835675434e2 * pow(x, 0.12e2) * pow(L, 0.28e2) + 0.8203174605e2 * pow(x, 0.8e1) * pow(L, 0.32e2) - 0.8685714297e2 * pow(x, 0.6e1) * pow(L, 0.34e2) + 0.5066666672e2 * pow(x, 0.4e1) * pow(L, 0.36e2) - 0.4772756133e2 * pow(x, 0.10e2) * pow(L, 0.30e2) - 0.4895134493e1 * pow(x, 0.14e2) * pow(L, 0.26e2) - 0.3278531612e-16 * pow(x, 0.38e2) * L * L + 0.6073479810e-14 * pow(x, 0.36e2) * pow(L, 0.4e1) - 0.6741562588e-12 * pow(x, 0.34e2) * pow(L, 0.6e1) + 0.5014037176e-10 * pow(x, 0.32e2) * pow(L, 0.8e1) - 0.2647411628e-8 * pow(x, 0.30e2) * pow(L, 0.10e2) + 0.1025872006e-6 * pow(x, 0.28e2) * pow(L, 0.12e2) - 0.2975028817e-5 * pow(x, 0.26e2) * pow(L, 0.14e2) + 0.7996418575e-19 * pow(x, 0.40e2));
  return retval;
}

double HO_42( const double L, const double x )
{
  double retval = -0.9849562077e25 * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.85e2 / 0.2e1) * (0.9999999997e0 * pow(L, 0.42e2) - 0.1599283714e-18 * pow(x, 0.42e2) - 0.4200000001e2 * x * x * pow(L, 0.40e2) + 0.2800000002e3 * pow(x, 0.4e1) * pow(L, 0.38e2) - 0.2108551674e-3 * pow(x, 0.26e2) * pow(L, 0.16e2) + 0.3807107186e-2 * pow(x, 0.24e2) * pow(L, 0.18e2) - 0.5253807919e-1 * pow(x, 0.22e2) * pow(L, 0.20e2) + 0.5516498312e0 * pow(x, 0.20e2) * pow(L, 0.22e2) - 0.4367227831e1 * pow(x, 0.18e2) * pow(L, 0.24e2) - 0.1101405260e3 * pow(x, 0.14e2) * pow(L, 0.28e2) - 0.6890666666e3 * pow(x, 0.10e2) * pow(L, 0.32e2) + 0.9119999997e3 * pow(x, 0.8e1) * pow(L, 0.34e2) - 0.7093333339e3 * pow(x, 0.6e1) * pow(L, 0.36e2) + 0.3340929289e3 * pow(x, 0.12e2) * pow(L, 0.30e2) + 0.2569945608e2 * pow(x, 0.16e2) * pow(L, 0.26e2) + 0.6884916383e-16 * pow(x, 0.40e2) * L * L - 0.1342558695e-13 * pow(x, 0.38e2) * pow(L, 0.4e1) + 0.1573031271e-11 * pow(x, 0.36e2) * pow(L, 0.6e1) - 0.1238762125e-9 * pow(x, 0.34e2) * pow(L, 0.8e1) + 0.6949455521e-8 * pow(x, 0.32e2) * pow(L, 0.10e2) - 0.2872441615e-6 * pow(x, 0.30e2) * pow(L, 0.12e2) + 0.8925086448e-5 * pow(x, 0.28e2) * pow(L, 0.14e2));
  return retval;
}

double HO_43( const double L, const double x )
{
  double retval = -0.5989635237e27 * x * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.87e2 / 0.2e1) * (0.9999999998e0 * pow(L, 0.42e2) - 0.3719264448e-20 * pow(x, 0.42e2) - 0.1400000001e2 * x * x * pow(L, 0.40e2) + 0.5600000002e2 * pow(x, 0.4e1) * pow(L, 0.38e2) - 0.7809450640e-5 * pow(x, 0.26e2) * pow(L, 0.16e2) + 0.1522842874e-3 * pow(x, 0.24e2) * pow(L, 0.18e2) - 0.2284264310e-2 * pow(x, 0.22e2) * pow(L, 0.20e2) + 0.2626903957e-1 * pow(x, 0.20e2) * pow(L, 0.22e2) - 0.2298540961e0 * pow(x, 0.18e2) * pow(L, 0.24e2) - 0.7342701731e1 * pow(x, 0.14e2) * pow(L, 0.28e2) - 0.6264242416e2 * pow(x, 0.10e2) * pow(L, 0.32e2) + 0.1013333333e3 * pow(x, 0.8e1) * pow(L, 0.34e2) - 0.1013333333e3 * pow(x, 0.6e1) * pow(L, 0.36e2) + 0.2569945605e2 * pow(x, 0.12e2) * pow(L, 0.30e2) + 0.1511732710e1 * pow(x, 0.16e2) * pow(L, 0.26e2) + 0.1679247897e-17 * pow(x, 0.40e2) * L * L - 0.3442458191e-15 * pow(x, 0.38e2) * pow(L, 0.4e1) + 0.4251435864e-13 * pow(x, 0.36e2) * pow(L, 0.6e1) - 0.3539320356e-11 * pow(x, 0.34e2) * pow(L, 0.8e1) + 0.2105895611e-9 * pow(x, 0.32e2) * pow(L, 0.10e2) - 0.9265940688e-8 * pow(x, 0.30e2) * pow(L, 0.12e2) + 0.3077616016e-6 * pow(x, 0.28e2) * pow(L, 0.14e2));
  return retval;
}

double HO_44( const double L, const double x )
{
  double retval = 0.4235311691e27 * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.89e2 / 0.2e1) * (-0.4400000004e2 * x * x * pow(L, 0.42e2) + 0.3080000002e3 * pow(x, 0.4e1) * pow(L, 0.40e2) - 0.8213333333e3 * pow(x, 0.6e1) * pow(L, 0.38e2) + 0.2454398773e-4 * pow(x, 0.28e2) * pow(L, 0.16e2) - 0.5154237421e-3 * pow(x, 0.26e2) * pow(L, 0.18e2) + 0.8375635805e-2 * pow(x, 0.24e2) * pow(L, 0.20e2) - 0.1050761583e0 * pow(x, 0.22e2) * pow(L, 0.22e2) + 0.1011358023e1 * pow(x, 0.20e2) * pow(L, 0.24e2) + 0.4038485952e2 * pow(x, 0.16e2) * pow(L, 0.28e2) + 0.4593777771e3 * pow(x, 0.12e2) * pow(L, 0.32e2) - 0.8917333325e3 * pow(x, 0.10e2) * pow(L, 0.34e2) + 0.1114666667e4 * pow(x, 0.8e1) * pow(L, 0.36e2) - 0.1615394380e3 * pow(x, 0.14e2) * pow(L, 0.30e2) - 0.7390693245e1 * pow(x, 0.18e2) * pow(L, 0.26e2) - 0.3518424167e-17 * pow(x, 0.42e2) * L * L + 0.7573408020e-15 * pow(x, 0.40e2) * pow(L, 0.4e1) - 0.9845430423e-13 * pow(x, 0.38e2) * pow(L, 0.6e1) + 0.8651671985e-11 * pow(x, 0.36e2) * pow(L, 0.8e1) - 0.5450553348e-9 * pow(x, 0.34e2) * pow(L, 0.10e2) + 0.2548133690e-7 * pow(x, 0.32e2) * pow(L, 0.12e2) - 0.9027673648e-6 * pow(x, 0.30e2) * pow(L, 0.14e2) + 0.9999999999e0 * pow(L, 0.44e2) + 0.7438528896e-20 * pow(x, 0.44e2));
  return retval;
}

double HO_45( const double L, const double x )
{
  double retval = 0.2695335856e29 * x * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.91e2 / 0.2e1) * (-0.1466666668e2 * x * x * pow(L, 0.42e2) + 0.6159999995e2 * pow(x, 0.4e1) * pow(L, 0.40e2) - 0.1173333333e3 * pow(x, 0.6e1) * pow(L, 0.38e2) + 0.8463444042e-6 * pow(x, 0.28e2) * pow(L, 0.16e2) - 0.1908976822e-4 * pow(x, 0.26e2) * pow(L, 0.18e2) + 0.3350254320e-3 * pow(x, 0.24e2) * pow(L, 0.20e2) - 0.4568528621e-2 * pow(x, 0.22e2) * pow(L, 0.22e2) + 0.4815990587e-1 * pow(x, 0.20e2) * pow(L, 0.24e2) + 0.2375579970e1 * pow(x, 0.16e2) * pow(L, 0.28e2) + 0.3533675205e2 * pow(x, 0.12e2) * pow(L, 0.32e2) - 0.8106666652e2 * pow(x, 0.10e2) * pow(L, 0.34e2) + 0.1238518517e3 * pow(x, 0.8e1) * pow(L, 0.36e2) - 0.1076929587e2 * pow(x, 0.14e2) * pow(L, 0.30e2) - 0.3889838547e0 * pow(x, 0.18e2) * pow(L, 0.26e2) - 0.8182381781e-19 * pow(x, 0.42e2) * L * L + 0.1847172687e-16 * pow(x, 0.40e2) * pow(L, 0.4e1) - 0.2524469338e-14 * pow(x, 0.38e2) * pow(L, 0.6e1) + 0.2338289725e-12 * pow(x, 0.36e2) * pow(L, 0.8e1) - 0.1557300956e-10 * pow(x, 0.34e2) * pow(L, 0.10e2) + 0.7721617239e-9 * pow(x, 0.32e2) * pow(L, 0.12e2) - 0.2912152787e-7 * pow(x, 0.30e2) * pow(L, 0.14e2) + 0.1000000000e1 * pow(L, 0.44e2) + 0.1653006420e-21 * pow(x, 0.44e2));
  return retval;
}

double HO_46( const double L, const double x )
{
  double retval = -0.1905890261e29 * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.93e2 / 0.2e1) * (0.9999999999e0 * pow(L, 0.46e2) + 0.3373333331e3 * pow(x, 0.4e1) * pow(L, 0.42e2) - 0.9445333326e3 * pow(x, 0.6e1) * pow(L, 0.40e2) + 0.1349333332e4 * pow(x, 0.8e1) * pow(L, 0.38e2) - 0.2595456171e-5 * pow(x, 0.30e2) * pow(L, 0.16e2) + 0.6272352414e-4 * pow(x, 0.28e2) * pow(L, 0.18e2) - 0.1185474606e-2 * pow(x, 0.26e2) * pow(L, 0.20e2) + 0.1751269304e-1 * pow(x, 0.24e2) * pow(L, 0.22e2) - 0.2013959700e0 * pow(x, 0.22e2) * pow(L, 0.24e2) - 0.1214185318e2 * pow(x, 0.18e2) * pow(L, 0.28e2) - 0.2322129422e3 * pow(x, 0.14e2) * pow(L, 0.32e2) + 0.6215111096e3 * pow(x, 0.12e2) * pow(L, 0.34e2) - 0.1139437035e4 * pow(x, 0.10e2) * pow(L, 0.36e2) + 0.6192345126e2 * pow(x, 0.16e2) * pow(L, 0.30e2) + 0.1789325733e1 * pow(x, 0.20e2) * pow(L, 0.26e2) + 0.1710861645e-18 * pow(x, 0.44e2) * L * L - 0.4046187790e-16 * pow(x, 0.42e2) * pow(L, 0.4e1) + 0.5806279478e-14 * pow(x, 0.40e2) * pow(L, 0.6e1) - 0.5661122490e-12 * pow(x, 0.38e2) * pow(L, 0.8e1) + 0.3979769109e-10 * pow(x, 0.36e2) * pow(L, 0.10e2) - 0.2089378782e-8 * pow(x, 0.34e2) * pow(L, 0.12e2) + 0.8372439261e-7 * pow(x, 0.32e2) * pow(L, 0.14e2) - 0.4600000004e2 * x * x * pow(L, 0.44e2) - 0.3306012840e-21 * pow(x, 0.46e2));
  return retval;
}

double HO_47( const double L, const double x )
{
  double retval = -0.1266807853e31 * x * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.95e2 / 0.2e1) * (0.1000000000e1 * pow(L, 0.46e2) + 0.6746666656e2 * pow(x, 0.4e1) * pow(L, 0.42e2) - 0.1349333332e3 * pow(x, 0.6e1) * pow(L, 0.40e2) + 0.1499259255e3 * pow(x, 0.8e1) * pow(L, 0.38e2) - 0.8372439257e-7 * pow(x, 0.30e2) * pow(L, 0.16e2) + 0.2162880141e-5 * pow(x, 0.28e2) * pow(L, 0.18e2) - 0.4390646686e-4 * pow(x, 0.26e2) * pow(L, 0.20e2) + 0.7005077210e-3 * pow(x, 0.24e2) * pow(L, 0.22e2) - 0.8756346516e-2 * pow(x, 0.22e2) * pow(L, 0.24e2) - 0.6390449039e0 * pow(x, 0.18e2) * pow(L, 0.28e2) - 0.1548086280e2 * pow(x, 0.14e2) * pow(L, 0.32e2) + 0.4780854687e2 * pow(x, 0.12e2) * pow(L, 0.34e2) - 0.1035851848e3 * pow(x, 0.10e2) * pow(L, 0.36e2) + 0.3642555951e1 * pow(x, 0.16e2) * pow(L, 0.30e2) + 0.8520598719e-1 * pow(x, 0.20e2) * pow(L, 0.26e2) + 0.3801914764e-20 * pow(x, 0.44e2) * L * L - 0.9409739040e-18 * pow(x, 0.42e2) * pow(L, 0.4e1) + 0.1416165725e-15 * pow(x, 0.40e2) * pow(L, 0.6e1) - 0.1451569868e-13 * pow(x, 0.38e2) * pow(L, 0.8e1) + 0.1075613272e-11 * pow(x, 0.36e2) * pow(L, 0.10e2) - 0.5969653658e-10 * pow(x, 0.34e2) * pow(L, 0.12e2) + 0.2537102805e-8 * pow(x, 0.32e2) * pow(L, 0.14e2) - 0.1533333330e2 * x * x * pow(L, 0.44e2) - 0.7034069867e-23 * pow(x, 0.46e2));
  return retval;
}

double HO_48( const double L, const double x )
{
  double retval = 0.8957684231e30 * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.97e2 / 0.2e1) * (-0.4799999989e2 * x * x * pow(L, 0.46e2) - 0.1079466666e4 * pow(x, 0.6e1) * pow(L, 0.42e2) + 0.1619199996e4 * pow(x, 0.8e1) * pow(L, 0.40e2) - 0.1439288885e4 * pow(x, 0.10e2) * pow(L, 0.38e2) + 0.2511731777e-6 * pow(x, 0.32e2) * pow(L, 0.16e2) - 0.6921216451e-5 * pow(x, 0.30e2) * pow(L, 0.18e2) + 0.1505364578e-3 * pow(x, 0.28e2) * pow(L, 0.20e2) - 0.2586490047e-2 * pow(x, 0.26e2) * pow(L, 0.22e2) + 0.3502538605e-1 * pow(x, 0.24e2) * pow(L, 0.24e2) + 0.3067415540e1 * pow(x, 0.20e2) * pow(L, 0.28e2) + 0.9288517677e2 * pow(x, 0.16e2) * pow(L, 0.32e2) - 0.3278300358e3 * pow(x, 0.14e2) * pow(L, 0.34e2) + 0.8286814789e3 * pow(x, 0.12e2) * pow(L, 0.36e2) - 0.1942696508e2 * pow(x, 0.18e2) * pow(L, 0.30e2) - 0.3718079443e0 * pow(x, 0.22e2) * pow(L, 0.26e2) - 0.7934430813e-20 * pow(x, 0.46e2) * L * L + 0.2053033973e-17 * pow(x, 0.44e2) * pow(L, 0.4e1) - 0.3236950229e-15 * pow(x, 0.42e2) * pow(L, 0.6e1) + 0.3483767685e-13 * pow(x, 0.40e2) * pow(L, 0.8e1) - 0.2717338793e-11 * pow(x, 0.38e2) * pow(L, 0.10e2) + 0.1591907643e-9 * pow(x, 0.36e2) * pow(L, 0.12e2) - 0.7163584391e-8 * pow(x, 0.34e2) * pow(L, 0.14e2) + 0.3679999994e3 * pow(x, 0.4e1) * pow(L, 0.44e2) + 0.1000000000e1 * pow(L, 0.48e2) + 0.1406813973e-22 * pow(x, 0.48e2));
  return retval;
}

double HO_49( const double L, const double x )
{
  double retval = 0.6207358463e32 * x * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.99e2 / 0.2e1) * (-0.1600000002e2 * x * x * pow(L, 0.46e2) - 0.1542095239e3 * pow(x, 0.6e1) * pow(L, 0.42e2) + 0.1799111111e3 * pow(x, 0.8e1) * pow(L, 0.40e2) - 0.1308444443e3 * pow(x, 0.10e2) * pow(L, 0.38e2) + 0.7611308433e-8 * pow(x, 0.32e2) * pow(L, 0.16e2) - 0.2232650474e-6 * pow(x, 0.30e2) * pow(L, 0.18e2) + 0.5190912350e-5 * pow(x, 0.28e2) * pow(L, 0.20e2) - 0.9579592789e-4 * pow(x, 0.26e2) * pow(L, 0.22e2) + 0.1401015445e-2 * pow(x, 0.24e2) * pow(L, 0.24e2) + 0.1460674070e0 * pow(x, 0.20e2) * pow(L, 0.28e2) + 0.5463833941e1 * pow(x, 0.16e2) * pow(L, 0.32e2) - 0.2185533577e2 * pow(x, 0.14e2) * pow(L, 0.34e2) + 0.6374472931e2 * pow(x, 0.12e2) * pow(L, 0.36e2) - 0.1022471849e1 * pow(x, 0.18e2) * pow(L, 0.30e2) - 0.1616556283e-1 * pow(x, 0.22e2) * pow(L, 0.26e2) - 0.1688176773e-21 * pow(x, 0.46e2) * L * L + 0.4562297728e-19 * pow(x, 0.44e2) * pow(L, 0.4e1) - 0.7527791248e-17 * pow(x, 0.42e2) * pow(L, 0.6e1) + 0.8496994373e-15 * pow(x, 0.40e2) * pow(L, 0.8e1) - 0.6967535385e-13 * pow(x, 0.38e2) * pow(L, 0.10e2) + 0.4302453099e-11 * pow(x, 0.36e2) * pow(L, 0.12e2) - 0.2046738402e-9 * pow(x, 0.34e2) * pow(L, 0.14e2) + 0.7360000011e2 * pow(x, 0.4e1) * pow(L, 0.44e2) + 0.1000000000e1 * pow(L, 0.48e2) + 0.2871048930e-24 * pow(x, 0.48e2));
  return retval;
}

double HO_50( const double L, const double x )
{
  double retval = -0.4389265261e32 * exp(-0.5000000000e0 * x * x * pow(L, -0.2e1)) * pow(L, -0.101e3 / 0.2e1) * (0.4000000009e3 * pow(x, 0.4e1) * pow(L, 0.46e2) + 0.1927619048e4 * pow(x, 0.8e1) * pow(L, 0.42e2) - 0.1799111110e4 * pow(x, 0.10e2) * pow(L, 0.40e2) + 0.1090370370e4 * pow(x, 0.12e2) * pow(L, 0.38e2) - 0.2238620128e-7 * pow(x, 0.34e2) * pow(L, 0.16e2) + 0.6977032732e-6 * pow(x, 0.32e2) * pow(L, 0.18e2) - 0.1730304117e-4 * pow(x, 0.30e2) * pow(L, 0.20e2) + 0.3421283139e-3 * pow(x, 0.28e2) * pow(L, 0.22e2) - 0.5388520945e-2 * pow(x, 0.26e2) * pow(L, 0.24e2) - 0.6639427594e0 * pow(x, 0.22e2) * pow(L, 0.28e2) - 0.3035463302e2 * pow(x, 0.18e2) * pow(L, 0.32e2) + 0.1365958486e3 * pow(x, 0.16e2) * pow(L, 0.34e2) - 0.4553194953e3 * pow(x, 0.14e2) * pow(L, 0.36e2) + 0.5112359247e1 * pow(x, 0.20e2) * pow(L, 0.30e2) + 0.6735651180e-1 * pow(x, 0.24e2) * pow(L, 0.26e2) + 0.3517034944e-21 * pow(x, 0.48e2) * L * L - 0.9918038541e-19 * pow(x, 0.46e2) * pow(L, 0.4e1) + 0.1710861648e-16 * pow(x, 0.44e2) * pow(L, 0.6e1) - 0.2023093899e-14 * pow(x, 0.42e2) * pow(L, 0.8e1) + 0.1741883847e-12 * pow(x, 0.40e2) * pow(L, 0.10e2) - 0.1132224500e-10 * pow(x, 0.38e2) * pow(L, 0.12e2) + 0.5685384450e-9 * pow(x, 0.36e2) * pow(L, 0.14e2) - 0.1226666668e4 * pow(x, 0.6e1) * pow(L, 0.44e2) - 0.5000000006e2 * x * x * pow(L, 0.48e2) - 0.5742097860e-24 * pow(x, 0.50e2) + 0.1000000000e1 * pow(L, 0.50e2));
  return retval;
}

TEF EF_HO[] = {&HO_0, &HO_1, &HO_2, &HO_3, &HO_4, &HO_5, &HO_6, &HO_7, &HO_8, &HO_9, &HO_10, &HO_11, &HO_12, &HO_13, &HO_14, &HO_15, &HO_16, &HO_17, &HO_18, &HO_19, &HO_20, &HO_21, &HO_22, &HO_23, &HO_24, &HO_25, &HO_26, &HO_27, &HO_28, &HO_29, &HO_30, &HO_31, &HO_32, &HO_33, &HO_34, &HO_35, &HO_36, &HO_37, &HO_38, &HO_39, &HO_40, &HO_41, &HO_42, &HO_43, &HO_44, &HO_45, &HO_46, &HO_47, &HO_48, &HO_49, &HO_50};
#endif
