%This program is used to test the bead-spring model motion;


rp = RandomWalkParams('dimension',3,'simulation',1);

rm = RouseModel('rp',rp,'b',1,'fixedBeads',[],'connectedBeads',[]);

rm.Simulation

rm.Plot