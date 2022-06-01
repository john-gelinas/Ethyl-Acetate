function rmin = rminCalc(aa,ab,ac,ad,xa,xb,xc,xd,splittype)
%for splittype 1-3: aa = aac, ab = abc, ac = acc = 1
%for splittype 4-6: aa = aad, ab = abd, ac = acd, ad = add = 1
%splittype 1: A/BC
%splittype 2: AB/BC^a
%splittype 3: AB/C
%splittype 4: A/BCD
%splittype 5: AB/CD
%splittype 6: ABC/D
%splittype 7: AB/BCD
%splittype 8: ABC/CD
%splittype 9: ABC/BCD
%splittype 10: A/B

if splittype == 1
    rmin = ab.*(xa+xb)./(xa.*(aa-ab)) + xc./(xa.*(aa-1));
elseif splittype == 2
    error('split type not supported')
elseif splittype == 3
    rmin = ((xb+xc)./(ab-1) + xa./(aa-1))./((xa+xb).*(1+xa.*xc));
elseif splittype == 4
    rmin = ab.*(xa+xb)./(xa.*(aa-ab)) + ac.*xc./(xa.*(aa-ac)) + xd./(aa-1);
elseif splittype == 5
    rmin = (ac.*xa./(aa-ac) + ac*(xb+xc)./(ab-ac)) ./ ((xa+xb).*(1+xa.*(xc+xd))) ...
        + xd.*(xa./(aa-1) + xb./(ab-1))./((xa+xb).^2);
elseif splittype == 6
    rmin = (xa./(aa-1) + xb./(ab-1) + (xc+xd)./(ac-1))./((1-xd).*(1+xd.*(xa+xb)));
elseif splittype == 10
    rmin = 1./((aa-1).*xa);
else
    error('split type not supported')
end
end
    
    
    