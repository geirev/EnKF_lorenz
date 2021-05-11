module m_fex
contains
subroutine fex(neq,t,y,ydot)
   use mod_lorenzpar
   implicit none
   integer neq
   real y(neq),ydot(neq),t

   ydot(1)=-sigma*y(1)+sigma*y(2)
   ydot(2)=-y(1)*y(3)+r*y(1)-y(2)
   ydot(3)=y(1)*y(2)-b*y(3)

end subroutine fex
end module m_fex
