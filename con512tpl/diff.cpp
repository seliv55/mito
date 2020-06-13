double o2deliv(double pm){
  double kvp=kv*pat,kv2pm=kv2*pm;
  double kv2=kv*kv, kt2=kt*kt, kl2=kl*kl,km2=km*km, kvpm=kv*pm,kp2=kp*kp,pm2=pm*pm,pat2=pat*pat,vo22=vo2a*vo2a;
 pal = (kvp - vo2a)/kv;
 pa = (-kvp + kl*kvp - kl*kvpm + kl*kt*kvpm + vo2a - kl*vo2a)/((-1 + kl*kt)*kv);
 pv = ((-kt)*kvp + kl*kt*kvp - kvpm + kt*kvpm + kt*vo2a - kl*kt*vo2a)/((-1 + kl*kt)*kv);
 vo2b = -((com*k0*km*kv*vb*(-kv*pat + kl*kvp + kt*kvp - kl2*kt*kvp - kl*kt2*kvp + kl2*kt2*kvp + kvpm - kl*kvpm - kt*kvpm + kl2*kt*kvpm + kl*kt2*kvpm - kl2*kt2*kvpm + vo2a - kl*vo2a - kt*vo2a + kl2*kt*vo2a + kl*kt2*vo2a - kl2*kt2*vo2a))/(km2*kv2 - 2*km*kp*kv2 + kp2*kv2 - 2*kl*km2*kt*kv2 + 4*kl*km*kp*kt*kv2 - 2*kl*kp2*kt*kv2 + kl2*km2*kt2*kv2 - 2*kl2*km*kp*kt2*kv2 + kl2*kp2*kt2*kv2 + km*kv2*pat - kl*km*kv2*pat - kp*kv2*pat + kl*kp*kv2*pat + km*kt*kv2*pat - 2*kl*km*kt*kv2*pat + kl2*km*kt*kv2*pat - kp*kt*kv2*pat + 2*kl*kp*kt*kv2*pat - kl2*kp*kt*kv2*pat - kl*km*kt2*kv2*pat + kl2*km*kt2*kv2*pat + kl*kp*kt2*kv2*pat - kl2*kp*kt2*kv2*pat + kt*kv2*pat2 - 2*kl*kt*kv2*pat2 + kl2*kt*kv2*pat2 + km*kv2pm + kl*km*kv2pm - kp*kv2pm - kl*kp*kv2pm - km*kt*kv2pm - 2*kl*km*kt*kv2pm - 
kl2*km*kt*kv2pm + kp*kt*kv2pm + 2*kl*kp*kt*kv2pm + kl2*kp*kt*kv2pm + kl*km*kt2*kv2pm + kl2*km*kt2*kv2pm - kl*kp*kt2*kv2pm - kl2*kp*kt2*kv2pm + kv2*pat*pm - kl*kv2*pat*pm - kt*kv2*pat*pm + 2*kl*kt*kv2*pat*pm - kl2*kt*kv2*pat*pm - kl*kt2*kv2*pat*pm + kl2*kt2*kv2*pat*pm + kl*kv2*pm2 - 2*kl*kt*kv2*pm2 + kl*kt2*kv2*pm2 - km*kv*vo2a + kl*km*kv*vo2a + kp*kv*vo2a - kl*kp*kv*vo2a - km*kt*kv*vo2a + 2*kl*km*kt*kv*vo2a - kl2*km*kt*kv*vo2a + kp*kt*kv*vo2a -        2*kl*kp*kt*kv*vo2a + kl2*kp*kt*kv*vo2a + kl*km*kt2*kv*vo2a - kl2*km*kt2*kv*vo2a - kl*kp*kt2*kv*vo2a + kl2*kp*kt2*kv*vo2a - 2*kt*kvp*vo2a + 4*kl*kt*kvp*vo2a - 2*kl2*kt*kvp*vo2a - kvpm*vo2a + kl*kvpm*vo2a + kt*kvpm*vo2a - 2*kl*kt*kvpm*vo2a + kl2*kt*kvpm*vo2a + kl*kt2*kvpm*vo2a - kl2*kt2*kvpm*vo2a + kt*vo22 - 2*kl*kt*vo22 + kl2*kt*vo22));
return vo2b;
}
