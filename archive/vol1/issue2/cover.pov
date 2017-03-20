// This is a scene file for POV-Ray 3.6 available at http://www.povray.org/
// This file was written by Edgar A. Bering IV for the Waterloo Mathematics Review,
// and is distributed under the same terms as the Review itsself

// Toggle high-quality rendering
#declare HQ=1;

#declare Bessel = function(n,s) { 
	select(s-.6,0.0,sqrt(2/(pi*s))*cos(s-pi/4-(pi*n)/2) ) 
	// The asymptotic approximation to the Bessel function
	// is good enough if you cut out the singularities at zero.
}
#declare n = 3; 
#declare jnk = 9.7610;  // j_3,2 from Abramowitz & Stegun


// Position the camera, set up focus and depth of field
camera {
	location <2.7*cos(pi/3),1 ,2.7*sin(pi/3)>
	look_at <0, 0, 0>
	right (8.5/11)*x
	up y
}

// Co-ordinate grid. Not drawn where the eigenfunction graph is to avoid ugliness
difference {
	plane {
		y,0
		pigment {
			brick rgb <1.0,1.0,1.0>,color transmit 1.0
			brick_size <.5,.5,.5> mortar .01
		}
		finish {
			diffuse 1
		}
	}
	union {
		difference {
			isosurface {
				function { Bessel(n,jnk*sqrt(pow(x,2)+pow(z,2)))*cos(n*atan2(z,x))-y }
				max_gradient 380
				contained_by { sphere { < 0,0,0>,1 } }
				#if(HQ>=2) accuracy .001 #end
				translate -.1*y
				pigment {
					//color rgb < 1.0,0,0>
					color transmit 1.0
				}
			}
			isosurface {
				function { Bessel(n,jnk*sqrt(pow(x,2)+pow(z,2)))*cos(n*atan2(z,x))-y }
				max_gradient 380
				#if(HQ>=2) accuracy .001 #end
				contained_by { sphere { < 0,0,0>,1 } }
				translate .1*y
				pigment {
					//color rgb < 0,1.0,0>
					color transmit 1.0
				}
			}
		}
		/*difference {
			isosurface {
				function { Bessel(n,jnk*sqrt(pow(x,2)+pow(z,2)))*cos(n*atan2(z,x))+y }
				max_gradient 380
				contained_by { sphere { < 0,0,0>,1 } }
				#if(HQ>=2) accuracy .001 #end
				translate .1*y
				pigment {
					//color rgb < 1.0,0,0>
					color transmit 1.0
				}
			}
			isosurface {
				function { Bessel(n,jnk*sqrt(pow(x,2)+pow(z,2)))*cos(n*atan2(z,x))+y }
				max_gradient 380
				#if(HQ>=2) accuracy .001 #end
				contained_by { sphere { < 0,0,0>,1 } }
				translate -.1*y
				pigment {
					//color rgb < 0,1.0,0>
					color transmit 1.0
				}
			}
		}*/
	}
}

// The eigenfunction
isosurface {
	function { Bessel(n,jnk*sqrt(pow(x,2)+pow(z,2)))*cos(n*atan2(z,x))-y }
	max_gradient 380
	#if(HQ>=2) accuracy .001 #end
	open
	contained_by { sphere { < 0,0,0>,1 } }
	pigment { function { (y-1)/2}
		color_map {
			[0.0, rgb <222/255,35/255,153/255>]
			[0.45, rgb <222/255,35/255,153/255>]
			//[0.33, rgb <1.0,93/255,197/255>]
			//[0.66, rgb <1.0,153/255,220/255>]
			[.8, rgb <1.0,1.0,1.0>]
			[1.0, rgb <1.0,1.0,1.0>]
		}
	}
	finish {
		ambient .25
	}
}

light_source { <-2, 4, 3> color rgb <1.0,1.0,1.0> spotlight radius 30 falloff 20 tightness 5 point_at <0,0,0> #if(HQ>=1) area_light <1,0,0>,<0,0,1>,6,6 adaptive 1 jitter  #end
}
light_source { <2, -4, -3> color rgb <1.0,1.0,1.0> }
