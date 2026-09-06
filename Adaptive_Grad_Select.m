function [GOSM,alpha,SM,FOSM] = Adaptive_Grad_Select(lambda,theta,GOSM,alpha,SM,FOSM)
lambda_sort=sort(lambda,'descend');

if lambda_sort(1)<=theta
    sum=lambda_sort(1);

    l=1;
   for j=1:size(lambda,1)-1

          if sum+lambda_sort(j+1)<= theta
              sum=sum+lambda_sort(j+1);
              l=l+1;
          else
              break;
          end
   end
     
      id=find( lambda > min( lambda_sort(1:l) )  );
      
      if ~isempty(find(id==1))
        GOSM=GOSM(:,id);
        alpha=alpha(:,id);
        SM=SM(:,id);
        FOSM=FOSM(id);
      else
             GOSM=GOSM(:,[1;id]);
             alpha=alpha(:,[1;id]);
             SM=SM(:,[1;id]);
             FOSM=FOSM([1;id]);
        
      end
   
else

              id=find( lambda ==  lambda_sort(1) ) ;

             if id==1

                 GOSM=GOSM(:,1);

                 SM=SM(:,1);

                 alpha=alpha(:, 1);

                 FOSM=FOSM(1);

             else

                 GOSM=GOSM(:,[1;id]);

                 SM=SM(:,[1;id]);

                 alpha=alpha(:, [1;id]);

                 FOSM=FOSM([1;id]);

             end







end




