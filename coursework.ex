defmodule CourseworkPi do
	@moduledoc false

	@default_eps 1.0e-6
	@recursion_limit 200_000

	def main(argv) do
		eps = parse_eps(argv)

		IO.puts("π approximation (Elixir)")
		IO.puts("eps=#{format_sci(eps)}")
		IO.puts("reference :math.pi() = #{fmt(:math.pi())}")
		IO.puts("")

		run_method("Leibniz", &leibniz_iter/1, &leibniz_rec/1, eps)
		run_method("Nilakantha", &nilakantha_iter/1, &nilakantha_rec/1, eps)
		run_method("Machin (arctan)", &machin_iter/1, &machin_rec/1, eps)

		:ok
	end

	defp run_method(name, iter_fun, rec_fun, eps) do
		{t_iter, {pi_i, n}} = :timer.tc(fn -> iter_fun.(eps) end)
		IO.puts("#{name} (iter): n=#{n}  pi=#{fmt(pi_i)}  err=#{fmt_err(pi_i)}  time=#{ms(t_iter)}ms")

		if n <= @recursion_limit do
			{t_rec, pi_r} = :timer.tc(fn -> rec_fun.(n) end)
			IO.puts("#{name} (rec):  n=#{n}  pi=#{fmt(pi_r)}  err=#{fmt_err(pi_r)}  time=#{ms(t_rec)}ms")
		else
			IO.puts("#{name} (rec):  skipped (n=#{n} > #{@recursion_limit})")
		end

		IO.puts("")
	end

	# =========================
	# Leibniz
	# π = 4 * Σ (-1)^k / (2k+1)
	# =========================

	def leibniz_iter(eps), do: series_iter(eps, 0.0, nil, 0, &leibniz_step/2)

	defp leibniz_step(k, sum) do
		term = (if rem(k, 2) == 0, do: 1.0, else: -1.0) / (2.0 * k + 1.0)
		sum2 = sum + term
		{4.0 * sum2, sum2}
	end

	def leibniz_rec(n) when n >= 1, do: 4.0 * leibniz_sum_rec(0, n - 1)
	defp leibniz_sum_rec(k, n) when k > n, do: 0.0
	defp leibniz_sum_rec(k, n), do: leibniz_term(k) + leibniz_sum_rec(k + 1, n)
	defp leibniz_term(k), do: (if rem(k, 2) == 0, do: 1.0, else: -1.0) / (2.0 * k + 1.0)

	# =========================
	# Nilakantha
	# π = 3 + Σ (-1)^(k+1) * 4 / ((2k)(2k+1)(2k+2)), k>=1
	# =========================

	def nilakantha_iter(eps), do: series_iter(eps, 0.0, nil, 0, &nilakantha_step/2)

	defp nilakantha_step(k0, sum) do
		k = k0 + 1
		a = 2.0 * k
		sign = if rem(k, 2) == 1, do: 1.0, else: -1.0
		term = sign * 4.0 / (a * (a + 1.0) * (a + 2.0))
		sum2 = sum + term
		{3.0 + sum2, sum2}
	end

	def nilakantha_rec(n) when n >= 1, do: 3.0 + nilakantha_sum_rec(1, n)
	defp nilakantha_sum_rec(k, n) when k > n, do: 0.0
	defp nilakantha_sum_rec(k, n), do: nilakantha_term(k) + nilakantha_sum_rec(k + 1, n)

	defp nilakantha_term(k) do
		a = 2.0 * k
		sign = if rem(k, 2) == 1, do: 1.0, else: -1.0
		sign * 4.0 / (a * (a + 1.0) * (a + 2.0))
	end

	# =========================
	# Machin formula
	# π = 16*arctan(1/5) - 4*arctan(1/239)
	# arctan(x) = Σ (-1)^k x^(2k+1)/(2k+1)
	# =========================

	def machin_iter(eps) do
		x1 = 1.0 / 5.0
		x2 = 1.0 / 239.0
		state = {0.0, x1, 1.0, 0.0, x2, 1.0}
		series_iter(eps, state, nil, 0, fn k, st -> machin_step(k, st, x1, x2) end)
	end

	defp machin_step(k, {s1, pow1, sign1, s2, pow2, sign2}, x1, x2) do
		denom = 2.0 * k + 1.0
		s1_2 = s1 + sign1 * (pow1 / denom)
		s2_2 = s2 + sign2 * (pow2 / denom)
		pow1_2 = pow1 * x1 * x1
		pow2_2 = pow2 * x2 * x2
		sign1_2 = -sign1
		sign2_2 = -sign2
		pi = 16.0 * s1_2 - 4.0 * s2_2
		{pi, {s1_2, pow1_2, sign1_2, s2_2, pow2_2, sign2_2}}
	end

	def machin_rec(n) when n >= 1 do
		x1 = 1.0 / 5.0
		x2 = 1.0 / 239.0
		16.0 * arctan_sum_rec(x1, 0, n - 1) - 4.0 * arctan_sum_rec(x2, 0, n - 1)
	end

	defp arctan_sum_rec(_x, k, n) when k > n, do: 0.0

	defp arctan_sum_rec(x, k, n) do
		sign = if rem(k, 2) == 0, do: 1.0, else: -1.0
		term = sign * :math.pow(x, 2 * k + 1) / (2.0 * k + 1.0)
		term + arctan_sum_rec(x, k + 1, n)
	end

	defp series_iter(eps, state, prev_pi, k, step_fun) do
		{pi, state2} = step_fun.(k, state)

		if prev_pi != nil and abs(pi - prev_pi) < eps do
			{pi, k + 1}
		else
			series_iter(eps, state2, pi, k + 1, step_fun)
		end
	end

	defp parse_eps([]), do: @default_eps
	defp parse_eps([eps_str | _]) do
		case Float.parse(eps_str) do
			{eps, _} when eps > 0.0 -> eps
			_ -> @default_eps
		end
	end

	defp fmt(x), do: :io_lib.format(~c"~.15f", [x]) |> IO.iodata_to_binary()
	defp fmt_err(pi), do: :io_lib.format(~c"~.3e", [abs(pi - :math.pi())]) |> IO.iodata_to_binary()
	defp ms(us), do: :io_lib.format(~c"~.3f", [us / 1000.0]) |> IO.iodata_to_binary()
	defp format_sci(x), do: :io_lib.format(~c"~.2e", [x]) |> IO.iodata_to_binary()
end

CourseworkPi.main(System.argv())

