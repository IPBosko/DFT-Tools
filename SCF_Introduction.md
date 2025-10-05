# Introduction to the Self-Consistent Field (SCF) Method

The main workhorse of DFT implementations is the SCF method. In this technical introduction, we will cover the general structure of this method and discuss the main focus of these tutorials.

## Overview

Again, the goal is to solve the equations of DFT to obtain the electron density of a molecular system. The electron density is required because it contains all the necessary information that can be used in various applications. 

The general outline of any SCF program is presented in the figure 

<p align="center">
  <img src="https://github.com/IPBosko/DFT-Tutorials/blob/main/ImageFolder/SCF_Outline.png" width="35%" height="35%" />
</p>

In this tutorial, we will cover all the green, orange, and blue parts of the diagram. However, since the green and blue parts are mostly the same for any SCF, most attention should be paid to the middle (orange) part, as it will contain all the differences between various approximations used in DFT. By the way, what are these "types of approximations" I am referring to all the time? The answer is in the following section.

## Density Functional Approximations (DFAs)

The DFT implementations are constructed in a way that requires approximating a part of the total energy functional. The most successful DFT implementation is the Kohn-Sham (KS-DFT) method. In KS-DFT, the goal is to approximate the exchange-correlation functional. The variety of those approximations can be presented as the so-called Jacob's ladder of KS-DFT. 

<p align="center">
  <img width="850" height="654" alt="Jacob's Ladder of KS-DFT" src="https://github.com/user-attachments/assets/fe77b440-23ce-4403-9b2b-357a23158cfd" />
</p>

In this diagram, the development of density functional approximations is given as a climb towards increasing accuracy. This is achieved via altering the structure of DFAs, including new ingredients, empirical parameters, etc. 

## Conclusion

We will climb this ladder step by step, covering all the necessary technical aspects. The first stop is the Hartree-Fock method, which is perfect for an easy and understandable introduction to the general structure of the SCF procedure.
