## Darwin and Ciocca's curve
bias = 0.1
dc_curve(Δf; a = 155.0 + bias, k = 0.057, s = 20.0) = a + k*Δf*exp(-(Δf^2)/(2*s^2))
mistuned_f = mistunings .* 620.0
plot!(mistuned_f, dc_curve.(mistuned_f .- 620.0),
      linestyle=:dash,
      linewidth=0.25,
      color=:gray,
      )

dc_f_points = [550.0, 570.0:10.0:670.0..., 690.0]
dc_f0_points = [154.9, 155.1, 154.75, 154.6, 154.5, 154.65, 155.4, 155.65, 155.75, 155.6, 155.4, 155.25, 155.20]
dc_bar_b = [154.75, 154.9, 154.5, 154.3, 154.25, 154.55, 155.2, 155.5, 155.6, 155.4, 155.25, 155.0, 155.0]
dc_bar_t = [155.1, 155.25, 154.95, 154.8, 154.65, 154.75, 155.0, 155.75, 156.0, 155.75, 155.55, 155.45, 155.45]
err_bottom = [abs(x-y) for (x,y) in (dc_bar_b, dc_f0_points)]
err_top = [abs(x-y) for (x,y) in (dc_bar_t, dc_f0_points)]
#=
err_bottom = [0.1 for k in 1:length(dc_f_points)]
err_top = [0.2 for k in 1:length(dc_f_points)]
=#

scatter!(dc_f_points, dc_f0_points,
         color=:gray,
         yerror = (err_bottom, err_top),
         markershape=:square,
         markersize=8,
         markerstrokewidth=0.5,
         markerstrokecolor=:black,
         markeralpha=0.3,
         markercolor=:black,
         )
