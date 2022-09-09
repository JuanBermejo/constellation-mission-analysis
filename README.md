# Constellation mission analysis

The aim of this study case is to evaluate the expansion of an existing constellation of satellites using the open source library of Astrodynamics Poliastro [[1]][poliastro]. In particular, it will focus on determining which orbits for the new satellites would perform better and how to deploy such satellites to those orbits using low-thrust maneuvers from a parking orbit.

## Description of the project

Satellite constellations have been used in the space sector for a long time but the number of operating constellations and its size has grown considerably in recent years, particularly with the deployment of new GNSS and communications constellations as GALILEO or Starlink. The design and operation of satellites, whether they are part of a constellation or not, requires software tools and algorithms to simulate their orbits. One of the most extended is the SGP4 propagator [[2]][vallado], based on a simplified perturbation model. This and other algorithms make use of two-line elements sets (TLEs) [[3]][kelso] to describe the orbit of each satellite, which are widely used and published for most active and inactive objetcs in orbit.

When studing new satellites, particuarly in the case of constellations, it is useful to simulate their potential orbit alongside existing satellites. To do so requires the generation of TLEs for the new satellites, which is not a straightforward process [[4]][rockwood].

The idea of this project is to study the extension of an existing constellation to improve its coverage of a certain region. To do so, the Astrodinamics library [Poliastro][poliastro] will be used. The study will be divided in two parts: **operational orbit design** and **deployment analysis**. For the first part, alternatives for the orbit of new satellites will be studied, assessing the improvement in the constellation performance. As previously said, to be able to simulate the potential satellites alongside those already deployed, it will be necessary to generate TLEs, for which the approach presented in [[4]][rockwood] will be followed. In the second part, the transference from the parking orbit to the operational orbit using low-thrust maneuvers will be studied.

As a closing point of the project, it will be proposed to contribute to Poliastro with those tools developed during the project which can be useful for others users (in particular, the TLE conversion).

## Preliminary timeline

The project will start by the end of september and it will last 4 months. Part of the following tasks will overlap over time.

- Getting familiar with Poliastro and the problem proposed (2 weeks)
- Prepare a proposal of how the problem will be solved (2 weeks)
- Implementation of the required code (4-6 weeks)
- Perform analysis (2-4 weeks)
- Prepare report (2-4 weeks)


## License

All the project will be developed under a [MIT][mit] license.

## References


[1] Cano Rodriguez, J. L., & Martínez Garrido, J. (2022). poliastro (Version v0.17.0) [Computer software]. https://github.com/poliastro/poliastro/
[2] David Vallado, Paul Crawford, Richard Hujsak, and T.S. Kelso. Revisiting spacetrack report #3. AIAA/AAS Astrodynamics Specialist Conference and Exhibit, 2012. doi: 10.2514/6.2006-6753. URL https://arc.aiaa.org/doi/abs/10.2514/6.2006-6753
[3] T S Kelso. Norad two-line element set format, Dec 2019. URL https://www.celestrak.com/NORAD/documentation/tle-fmt.php
[4] ROCKWOOD, Troy; STEEGER, Greg; STEIN, Matthew. Generating Realistic Two-Line Elements for Notional Space Vehicles and Constellations. arXiv preprint arXiv:2203.04204, 2022. https://doi.org/10.48550/arXiv.2203.04204

[//]: # (These are reference links used in the text)

   [poliastro]: <https://github.com/poliastro/poliastro>
   [rockwood]: <https://doi.org/10.48550/arXiv.2203.04204>
   [mit]: <https://github.com/git/git-scm.com/blob/main/MIT-LICENSE.txt>
   [kelso]: <https://www.celestrak.com/NORAD/documentation/tle-fmt.php>
   [vallado]: < https://arc.aiaa.org/doi/abs/10.2514/6.2006-6753>