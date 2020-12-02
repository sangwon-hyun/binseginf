context("Test hit and run, Radial")

## .hit_run_next_point_radial is correct

test_that(".hit_run_next_point_radial works", {
  set.seed(15)
  y <- rnorm(10)
  obj <- binseginf::bsfs(y, 2)
  poly <- binseginf::polyhedra(obj)
  segments <- .segments(length(y), jump_vec = binseginf::jumps(obj))

  null_mat <- .compute_nullspace(segments)
  res <- .hit_run_next_point_radial(y, null_mat, poly)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == length(y))
})

test_that(".hit_run_next_point_radial gives a sample with the correct properties", {
  set.seed(10)
  y <- rnorm(10)
  obj <- binseginf::bsfs(y, 2)
  poly <- binseginf::polyhedra(obj)
  segments <- .segments(length(y), jump_vec = binseginf::jumps(obj))

  null_mat <- .compute_nullspace(segments)
  res <- .hit_run_next_point_radial(y, null_mat, poly)

  expect_true(sum(abs(.segment_means(y, segments) -
                        .segment_means(res, segments))) < 1e-6)
  expect_true(abs(.l2norm(y) - .l2norm(res)) < 1e-6)
})

###########################

## .sampler_hit_run_radial is correct

test_that(".sampler_hit_run_radial works", {
  set.seed(50)
  y <- rnorm(10)
  obj <- binseginf::bsfs(y, 2)
  poly <- binseginf::polyhedra(obj)
  segments <- .segments(length(y), jump_vec = binseginf::jumps(obj))

  res <- .sampler_hit_run_radial(y, segments, poly, num_samp = 25, burn_in = 10)

  expect_true(is.numeric(res))
  expect_true(is.matrix(res))
  expect_true(all(dim(res) == c(length(y), 25)))
})

test_that(".sampler_hit_run_radial preserves segment means", {
  set.seed(10)
  y <- rnorm(10)
  obj <- binseginf::bsfs(y, 2)
  poly <- binseginf::polyhedra(obj)
  segments <- .segments(length(y), jump_vec = binseginf::jumps(obj))
  seg_mean <- .segment_means(y, segments)

  res <- .sampler_hit_run_radial(y, segments, poly, num_samp = 25, burn_in = 10)

  bool_vec <- sapply(1:25, function(x){
    seg_mean2 <- .segment_means(res[,x], segments)
    sum(abs(seg_mean - seg_mean)) < 1e-6
  })

  expect_true(all(bool_vec))
})

test_that(".sampler_hit_run_radial preserves l2norm", {
  set.seed(30)
  y <- rnorm(10)
  obj <- binseginf::bsfs(y, 2)
  poly <- binseginf::polyhedra(obj)
  segments <- .segments(length(y), jump_vec = binseginf::jumps(obj))
  norm1 <- .l2norm(y)

  res <- .sampler_hit_run_radial(y, segments, poly, num_samp = 25, burn_in = 10)

  bool_vec <- sapply(1:25, function(x){
    norm2 <- .l2norm(res[,x])
    abs(norm1 - norm2) < 1e-6
  })

  expect_true(all(bool_vec))
})

test_that(".sampler_hit_run_radial are all in polyhedra", {
  set.seed(10)
  y <- rnorm(10)
  obj <- binseginf::bsfs(y, 2)
  poly <- binseginf::polyhedra(obj)
  segments <- .segments(length(y), jump_vec = binseginf::jumps(obj))
  seg_mean <- .segment_means(y, segments)

  res <- .sampler_hit_run_radial(y, segments, poly, num_samp = 25, burn_in = 10)

  bool_vec <- sapply(1:25, function(x){
    .try_polyhedra(res[,x], poly)
  })

  expect_true(all(bool_vec))
})
