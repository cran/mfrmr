sampling_vcov_fixture <- function(nested=FALSE,one=FALSE,incomplete=TRUE) {
  x<-mfrmr_example_gtheory_crossed
  x<-x[x$Person%in%head(unique(x$Person),8) & x$Rater%in%head(unique(x$Rater),3) &
    x$Task%in%head(unique(x$Task),4),]
  if(one)x<-aggregate(x[c('Content','Organization')],x[c('Person','Task')],mean)
  if(incomplete)x<-x[-c(1,8,23),]
  mfrm_g_study(x,c('Content','Organization'),method='minque0',
    rater=if(one)NULL else 'Rater',nesting=if(nested)c(Rater='Task')else NULL)
}

test_that("component sampling covariance matches independent Gaussian product moments", {
  for(type in c('one','crossed','nested')) {
    g<-sampling_vcov_fixture(type=='nested',type=='one')
    x<-g$data;N<-nrow(x);H<-diag(N)-1/N
    groups<-if(type=='one')list(x$Person,x$Task,seq_len(N))else if(type=='nested')
      list(x$Person,x$Task,interaction(x$Rater,x$Task),interaction(x$Person,x$Task),seq_len(N))else
        list(x$Person,x$Rater,x$Task,interaction(x$Person,x$Rater),
          interaction(x$Person,x$Task),interaction(x$Rater,x$Task),seq_len(N))
    K<-lapply(groups,function(g)outer(g,g,`==`)*1);k<-length(K)
    B<-lapply(K,function(x)H%*%x%*%H)
    inverse<-solve(crossprod(vapply(B,as.vector,numeric(N*N))))
    A<-lapply(seq_len(k),function(s)Reduce(`+`,Map(`*`,B,inverse[s,])))
    P<-matrix(0,k*k,k*k)
    for(u in seq_len(k))for(v in seq_len(k)) {
      AKu<-lapply(A,function(a)a%*%K[[u]])
      AKv<-lapply(A,function(a)a%*%K[[v]])
      P[,u+(v-1)*k]<-as.vector(outer(seq_len(k),seq_len(k),Vectorize(function(s,t)
        2*sum(AKu[[s]]*t(AKv[[t]])))))
    }
    q<-vapply(g$components,function(m)m[1,1],numeric(1))
    expected<-matrix(P%*%solve(diag(k*k)+P,as.vector(tcrossprod(q))),k)
    got<-vcov(g,score='Content',assumption='normal')
    expect_equal(as.vector(got),as.vector(expected),tolerance=1e-9)
    expect_equal(attr(got,'summary')$SE^2,diag(expected),tolerance=1e-9)
    # Validate unbiasedness at a population covariance, without random draws.
    truth<-seq_len(k)/k+.5
    Eqq<-(diag(k*k)+P)%*%as.vector(tcrossprod(truth))
    V<-Reduce(`+`,Map(`*`,K,truth));AV<-lapply(A,function(a)a%*%V)
    C<-outer(seq_len(k),seq_len(k),Vectorize(function(s,t)2*sum(AV[[s]]*t(AV[[t]]))))
    expect_equal(as.vector(P%*%solve(diag(k*k)+P,Eqq)),as.vector(C),tolerance=1e-9)
  }
})

test_that("balanced MINQUE and ANOVA sampling SEs match existing independent-MS formulas", {
  for(type in c('one','crossed','nested')) {
    g<-sampling_vcov_fixture(type=='nested',type=='one',FALSE)
    a<-mfrm_g_study(g$data,g$design$scores,method='anova',
      rater=if(type=='one')NULL else 'Rater',nesting=if(type=='nested')c(Rater='Task')else NULL)
    for(w in list(NULL,c(Content=.4,Organization=-.8))) {
      x<-vcov(g,weights=w,assumption='normal')
      y<-vcov(a,weights=w,assumption='normal')
      ci<-confint(a,weights=w,assumption='normal')
      expect_equal(as.vector(x),as.vector(y),tolerance=1e-9)
      expect_equal(attr(y,'summary')$SE,ci$SE,tolerance=1e-9)
      expect_equal(attr(y,'summary')$Estimate,ci$Estimate,tolerance=1e-9)
    }
    s<-c('Content','Organization')
    x<-vcov(g,score=s,assumption='normal');y<-vcov(a,score=s,assumption='normal')
    ci<-confint(a,score=s,method='regression_region',assumption='normal')
    expect_equal(as.vector(x),as.vector(y),tolerance=1e-9)
    expect_equal(attr(y,'summary')$SE,ci$SE,tolerance=1e-9)
    a$data<-NULL
    expect_identical(vcov(a,score=s,assumption='normal'),y)
  }
})

test_that("between-score moment correction uses coupled equations and preserves units", {
  g<-sampling_vcov_fixture()
  s<-c('Content','Organization')
  x<-vcov(g,score=s,assumption='normal')
  reversed<-vcov(g,score=rev(s),assumption='normal')
  expect_equal(as.vector(x),as.vector(reversed),tolerance=1e-10)
  # Independent 3k-dimensional quadratic-form covariance identifies all cross terms.
  groups<-.mfrm_mvgt_moment_groups(g);N<-nrow(g$data);k<-length(groups);H<-diag(N)-1/N
  K<-lapply(groups,function(g)outer(g,g,`==`)*1)
  B<-lapply(K,function(x)H%*%x%*%H)
  inverse<-solve(crossprod(vapply(B,as.vector,numeric(N*N))))
  A<-lapply(seq_len(k),function(s)Reduce(`+`,Map(`*`,B,inverse[s,])))
  theta1<-seq_len(k)/k+.5;theta2<-rev(theta1)*1.4;cov<-sqrt(theta1*theta2)*seq(-.5,.5,length.out=k)
  V1<-Reduce(`+`,Map(`*`,K,theta1));V2<-Reduce(`+`,Map(`*`,K,theta2));Vc<-Reduce(`+`,Map(`*`,K,cov))
  V<-rbind(cbind(V1,Vc),cbind(Vc,V2));Z<-matrix(0,N,N)
  QV<-lapply(A,function(a)rbind(cbind(Z,a/2),cbind(a/2,Z))%*%V)
  expected<-outer(seq_len(k),seq_len(k),Vectorize(function(s,t)2*sum(QV[[s]]*t(QV[[t]]))))
  P<-.mfrm_mvgt_moment_polynomial(groups,g$estimation$kernel_gram)
  m<-.mfrm_mvgt_moment_map(P,TRUE);sym<-m$symmetric
  truth<-c(((outer(theta1,theta2)+outer(theta2,theta1))/2)[sym$index],tcrossprod(cov)[sym$index])
  expect_equal(as.vector(sym$expand%*%(m$weights%*%m$moment%*%truth)),as.vector(expected),tolerance=1e-8)
  data<-g$data;data$Content<-data$Content*(-2);data$Organization<-data$Organization*3
  scaled<-vcov(mfrm_g_study(data,s,method='minque0'),score=s,assumption='normal')
  expect_equal(as.vector(scaled),as.vector(x)*36,tolerance=1e-9)
  expect_equal(attr(scaled,'summary')$Estimate,-6*attr(x,'summary')$Estimate,tolerance=1e-9)
  expect_equal(attr(scaled,'summary')$SE,6*attr(x,'summary')$SE,tolerance=1e-9)
})

test_that("sampling covariance retains nonidentifiability, negative variance and saved state", {
  g<-sampling_vcov_fixture(incomplete=FALSE)
  data<-g$data[g$data$Rater%in%head(unique(g$data$Rater),2),]
  a<-mfrm_g_study(data,g$design$scores)
  x<-vcov(a,score=g$design$scores,assumption='normal');tab<-attr(x,'summary')
  expect_true(is.na(tab$SE[tab$Source=='Rater']))
  expect_match(tab$SEStatus[tab$Source=='Rater'],'not identifiable')
  expect_true(all(is.finite(tab$SE[tab$Source!='Rater'])))
  ci<-confint(a,score=g$design$scores,method='regression_region',assumption='normal')
  expect_equal(tab$SE,ci$SE,tolerance=1e-9)
  before<-serialize(a,NULL)
  invisible(vcov(a,assumption='normal'))
  expect_identical(serialize(a,NULL),before)
  path<-tempfile(fileext='.rds');on.exit(unlink(path),add=TRUE)
  saveRDS(x,path);expect_identical(readRDS(path),x)
  # Status propagation is separate from the independently checked moment algebra.
  local_mocked_bindings(.mfrm_mvgt_moment_correct=function(map,...) {
    k<-map$symmetric$k;V<-diag(k);V[1,1]<--.1
    list(covariance=V,estimable=matrix(TRUE,k,k),se=c(NA,rep(1,k-1)))
  })
  negative<-vcov(sampling_vcov_fixture(),assumption='normal');tab<-attr(negative,'summary')
  expect_lt(tab$SamplingVariance[1],0);expect_true(is.na(tab$SE[1]))
  expect_match(tab$SEStatus[1],'Negative')
  expect_false(attr(negative,'diagnostics')$positive_semidefinite)
})

test_that("sampling covariance rejects incompatible sources and preserves named composites", {
  g<-sampling_vcov_fixture();w<-c(Content=.5,Organization=-.5)
  x<-vcov(g,weights=w,assumption='normal')
  y<-vcov(g,weights=-3*w,assumption='normal')
  expect_equal(as.vector(y),81*as.vector(x),tolerance=1e-9)
  expect_equal(attr(y,'summary')$SE,9*attr(x,'summary')$SE,tolerance=1e-9)
  expect_error(vcov(g),'normal')
  expect_error(vcov(g,score=g$design$scores,weights=w,assumption='normal'),'without weights')
  expect_error(vcov(g,score=c('Content','Content'),assumption='normal'),'distinct')
  expect_error(vcov(g,assumption='normal',unused=1),'must be empty')
  changed<-g;changed$design$counts[1]<-100
  expect_error(vcov(changed,assumption='normal'),'counts')
  changed<-g;changed$estimation$kernel_gram[1,1]<-0
  expect_error(vcov(changed,assumption='normal'),'equations')
  changed<-g;changed$components[[1]][1,1]<-1
  expect_error(vcov(changed,assumption='normal'),'quadratic products')
  changed<-g;changed$data<-NULL
  expect_error(vcov(changed,assumption='normal'),'assignments')
  reordered<-mfrm_g_study(g$data[nrow(g$data):1,],g$design$scores,method='minque0')
  expect_equal(as.vector(vcov(reordered,weights=w,assumption='normal')),as.vector(x),tolerance=1e-9)
  renamed<-mfrm_g_study(g$data,g$design$scores,method='minque0',
    facets=c(Assessor='Rater',Occasion='Task'))
  v<-vcov(renamed,weights=w,assumption='normal')
  expect_equal(as.vector(v),as.vector(x),tolerance=1e-9)
  expect_match(rownames(v)[2],'Assessor')
})

test_that("unequal nested assignments and changed score units retain the moment target", {
  g<-sampling_vcov_fixture(TRUE,FALSE,FALSE);data<-g$data
  data<-data[!(data$Task==data$Task[1] & data$Rater==tail(unique(data$Rater),1)),]
  g<-mfrm_g_study(data,g$design$scores,method='minque0',nesting=c(Rater='Task'))
  expect_true(g$design$complete);expect_false(g$design$balanced)
  x<-vcov(g,score=g$design$scores,assumption='normal')
  expect_true(all(is.finite(attr(x,'summary')$SE)))
  swapped<-mfrm_g_study(data,g$design$scores,method='minque0',
    facets=c(Task='Task',Rater='Rater'),nesting=c(Rater='Task'))
  expect_equal(as.vector(vcov(swapped,score=g$design$scores,assumption='normal')),
    as.vector(x),tolerance=1e-8)
  data$Content<-data$Content*1e25;data$Organization<-data$Organization*1e-25
  scaled<-mfrm_g_study(data,g$design$scores,method='minque0',nesting=c(Rater='Task'))
  expect_equal(as.vector(vcov(scaled,score=g$design$scores,assumption='normal')),
    as.vector(x),tolerance=1e-8)
  zero<-g$data;zero$Content<-0
  constant<-mfrm_g_study(zero,g$design$scores,method='minque0',nesting=c(Rater='Task'))
  ans<-vcov(constant,score='Content',assumption='normal')
  expect_true(all(attr(ans,'summary')$SE==0))
  expect_true(all(attr(ans,'summary')$SEStatus=='Zero sampling-variance estimate'))
})
