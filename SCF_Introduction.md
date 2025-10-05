# Introduction to the Self-Consistent Field (SCF) Method

The main workhorse of DFT implementations is the SCF method. In this technical introduction, we will cover the general structure of this method and discuss the main focus of these tutorials.

## Overview

Again, the goal is to solve the equations of DFT to obtain the electron density of a molecular system. The electron density is required because it contains all the necessary information that can be used in various applications. 

The general outline of any SCF program is presented in the figure 

<p align="center">
  <img src="https://github.com/IPBosko/DFT-Tutorials/blob/main/ImageFolder/SCF_Outline.png" width="35%" height="35%" />
</p>

In this tutorial, we will cover all the green, orange, and blue parts of the diagram. However, since the green and blue parts are mostly the same for any SCF, most attention should be paid to the middle (orange) part, as it will contain all the differences between various approximations used in DFT. By the way, what are these "types of approximations" I am referring to all the time?

## Density Functional Approximations

