@echo off
cls
REM Parâmetros fixos
SET QtdeMaxGeracoes=500
SET TamanhoTorneio=2
SET QtdePontosCorte=2
SET parmLambda=0.01
SET RhoMinimo=0.2

setlocal EnableDelayedExpansion
REM G: quantidade de genes
for /l %%G in (10, 5, 20) do (	
	REM I: quantidade de indivíduos por geração
	for /l %%I in (10, 5, 20) do (	
		REM C: probabilidade de crossover
		for /l %%C in (70, 10, 100) do (	
			REM M: probabilidade de mutação
			for /l %%M in (5, 10, 30) do (	
			REM I: intensidade de mutação. importante: se alterar o 10 abaixo, deve alterar o código e recompilar
				for /l %%I in (1, 1, 10) do (
					Serial_novo %%G %QtdeMaxGeracoes% %%I %TamanhoTorneio% %%C %QtdePontosCorte% %%M %%I %parmLambda% %RhoMinimo% > Serial_Novo_parm_%%G_%%I_%%C_%%M_%%I.txt
				)
			)
		)
	)
)

pause