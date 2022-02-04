# RwR

RwR~ - A module for spatialization of audio source Released under GNU's GPL License.

Based on a concept by Shahrokh Yadegari in "space". Deformable polygon inner room by Yuki Minakawa and Julian Villegas.

RwR~ is a real-time implementation of a general model for the spatialization of audio sources. The central concept of the algorithm is that of a room within a room.

In this implementation, the outer room is assumed to be square, and the inner one is assumed to be a polygon. If a virtual source is in the listening room, the audio will be muted.

The inner room is the space delimited by the speakers which contain the listeners. The model simulates the behavior of the sound source within a user-defined virtual outer room, as heard from the inner room. The speakers act as "windows" through which sound from the outer room passes.

You can change the number of the windows from an argument of RwR~. To form a polygon, the number of the windows must be three or more.

inlets
----------
1 - audio source

2 - list of the windows in radius and angles

3 - the energy absorbed by the wall in the reflection in decibels

4 - X coordinate in meters of the audio source

5 - Y coordinate in meters of the audio source

6 - outer room size in meters
