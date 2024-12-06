# Programa simples para o entendimento do Movimento Harmônico Simples e Amortecido de uma partícula
# Destinado a alunos e tutores das disciplinas de Física Básica na Tutoria
# Não utilizamos métodos de integração como Runge-Kutta ou Preditor-Corretor, apenas a função analítica solução da EDO
# Desenvolvido por Msc.Antonio Kaeliton (antonio.kaeliton@estudante.ufjf.br) e Msc.Henrique Santiago (henrique.assis@estudante.ufjf.br)
# Disciplina: Tutoria - Doutorado em Física
# AVISO IMPORTANTÍSSIMO:
# CERTIFIQUE-SE DE COPIAR O PROGRAMA CORRETAMENTE, PARA EVITAR ERROS DE EXECUÇÃO
# NÃO ESQUEÇA DE CONFERIR O ARQUIVO DE ENTRADA PARA INSERIR OS VALORES DAS PRINCIPAIS CONSTANTES!
# LEIA O MANUAL PARA MAIS DÚVIDAS!

#Importando as principais bibliotecas

import math
import numpy as np
import matplotlib.pyplot as plt
#import plotly.graph_objects as go

def mhs_sa(xo, vo, m, k, npt, dt, wo, beta): 
    
    try:
        tempo = []
        posicao = []
        velocidade = []
        ener_cin = []
        ener_pot = []
        energia = []

        for t in np.arange(0, npt + 0.1, dt):

            x = math.exp(-beta * t) * (xo * math.cos(wo * t) + (vo/wo + beta * xo) * math.sin(wo * t))

            v = -beta * x + math.exp(-beta * t) * (-xo * wo * math.sin(wo * t) + (vo + beta * xo) * math.cos(wo * t))

            E_cin = 0.5 * m * v**2

            E_pot = 0.5 * k * x**2

            E_tot = E_cin + E_pot

            tempo.append(t)
            posicao.append(x)
            velocidade.append(v)
            ener_cin.append(E_cin)
            ener_pot.append(E_pot)
            energia.append(E_tot)
        
        return tempo, posicao, velocidade, ener_cin, ener_pot, energia

    except (ValueError, ZeroDivisionError, TypeError):
        print("Erro: Verifique seus dados de entrada")
        exit()

def mhs_sub(xo, vo, m, k, npt, dt, w, beta): 

    try:
        tempo = []
        posicao = []
        velocidade = []
        ener_cin = []
        ener_pot = []
        energia = []

        for t in np.arange(0, npt + 0.1, dt):

            x = math.exp(-beta * t) * (xo * math.cos(w * t) + (vo/w + beta * xo) * math.sin(w * t))

            v = -beta * x + math.exp(-beta * t) * (-xo * w * math.sin(w * t) + (vo + beta * xo) * math.cos(w * t))

            E_cin = 0.5 * m * v**2

            E_pot = 0.5 * k * x**2

            E_tot = E_cin + E_pot

            tempo.append(t)
            posicao.append(x)
            velocidade.append(v)
            ener_cin.append(E_cin)
            ener_pot.append(E_pot)
            energia.append(E_tot)
        
        return tempo, posicao, velocidade, ener_cin, ener_pot, energia
    
    except (ValueError, ZeroDivisionError, TypeError):
        print("Erro: Verifique seus dados de entrada")
        exit()

def mhs_cri(xo, vo, m, k, npt, dt, wo, beta): 

    try:
        tempo = []
        posicao = []
        velocidade = []
        ener_cin = []
        ener_pot = []
        energia = []

        for t in np.arange(0, npt + 0.1, dt):

            x = math.exp(-beta * t) * (xo + t/wo * (vo + beta * xo))
                        
            v = -beta * x + math.exp(-beta * t) * 1/wo * (vo + beta * xo)

            E_cin = 0.5 * m * v**2

            E_pot = 0.5 * k * x**2

            E_tot = E_cin + E_pot

            tempo.append(t)
            posicao.append(x)
            velocidade.append(v)
            ener_cin.append(E_cin)
            ener_pot.append(E_pot)
            energia.append(E_tot)
        
        return tempo, posicao, velocidade, ener_cin, ener_pot, energia
    
    except (ValueError, ZeroDivisionError, TypeError):
        print("Erro: Verifique seus dados de entrada")
        exit()

def mhs_sup(xo, vo, m, k, npt, dt, wo, beta): 

    try:
        alpha1 = beta + math.sqrt(beta**2 - wo**2)
        alpha2 = beta - math.sqrt(beta**2 - wo**2)

        C = (vo + alpha2 * xo) / (alpha2 - alpha1)
        D = xo - C

        tempo = []
        posicao = []
        velocidade = []
        ener_cin = []
        ener_pot = []
        energia = []
        
        for t in np.arange(0, npt + 0.1, dt):

            x = C * math.exp(-alpha1 * t) + D * math.exp(-alpha2 * t)
            
            v = -C * alpha1 * math.exp(-alpha1 * t) - D * alpha2 * math.exp(-alpha2 * t)
            
            E_cin = 0.5 * m * v**2

            E_pot = 0.5 * k * x**2

            E_tot = E_cin + E_pot

            tempo.append(t)
            posicao.append(x)
            velocidade.append(v)
            ener_cin.append(E_cin)
            ener_pot.append(E_pot)
            energia.append(E_tot)
        
        return tempo, posicao, velocidade, ener_cin, ener_pot, energia
    
    except (ValueError, ZeroDivisionError, TypeError):
        print("Erro: Verifique seus dados de entrada")
        exit()

try:
    with open('input_mhs.dat', 'r') as in_mhs:

        xo = float(in_mhs.readline().strip())
        vo = float(in_mhs.readline().strip())
        m = float(in_mhs.readline().strip())
        k = float(in_mhs.readline().strip())
        npt = float(in_mhs.readline().strip())
        dt = float(in_mhs.readline().strip())
        all = int(in_mhs.readline().strip())

    wo = math.sqrt(k/m)

except (ValueError, ZeroDivisionError):
    
    print("Erro: Verifique seus dados de entrada")
    exit()

if all == 0:

    try:

        mhs = int(input('Escolha o caso a ser avaliado: \n'
            '[0] Sem amortecimento; (Bônus -->) [1] Subamortecido; [2] Crítico; [3] Superamortecido \n            '))

        if (mhs not in (0 ,1 ,2 ,3)):
            print('Valor incorreto. Encerrando...')
            exit()

    except (ValueError, KeyboardInterrupt):
        print('Entrada inválida. Encerrando...') 
        exit()

    try:

        gr = int(input('Escolha o gráfico a ser plotado: \n'
            '[0] Posição x Tempo; [1] Velocidade x Tempo; [2] Energia x Tempo \n            '))

        if (gr not in (0 ,1 ,2)):
            print('Valor incorreto. Encerrando...')
            exit()
            
    except (ValueError, KeyboardInterrupt):
        print('Entrada inválida. Encerrando...') 
        exit()

######################################### Mude o valor de beta a partir daqui!. Atenção para as condições de cada caso
    if(mhs == 0):
    
        beta = 0.0
        w = wo
        print(f'Constante de amortecimento - Beta = {beta:.2f}')

        tempo, posicao, velocidade, ener_cin, ener_pot, energia = mhs_sa(xo, vo, m, k, npt, dt, wo, beta)

    elif(mhs == 1):
    
        beta = wo - 0.9 * wo 
        print(f'Constante de amortecimento - Beta = {beta:.2f}')
        
        if beta >= wo:
             print(f"Caso de Subamortecido - Erro: O valor de Beta ({beta:.2f}) deve ser menor que wo ({beta:.2f}).\n"
                   "Verifique seus Dados de Entrada")
             exit()

        w = math.sqrt(wo**2 - beta**2)
        tempo, posicao, velocidade, ener_cin, ener_pot, energia = mhs_sub(xo, vo, m, k, npt, dt, w, beta)

    elif(mhs == 2):
    
        beta = wo
        print(f'Constante de amortecimento - Beta = Omega = {beta:.2f}')

        tempo, posicao, velocidade, ener_cin, ener_pot, energia = mhs_cri(xo, vo, m, k, npt, dt, wo, beta)

    elif(mhs == 3):

        beta = wo + 0.9 * wo   
        print(f'Constante de amortecimento - Beta = {beta:.3f}')
        if beta <= wo:
             print(f"Caso de Superamortecido - Erro: o valor de Beta ({beta:.2f}) deve ser maior que wo ({wo:.2f}).\n"
                   "Verifique seus Dados de Entrada")
             exit()

        tempo, posicao, velocidade, ener_cin, ener_pot, energia = mhs_sup(xo, vo, m, k, npt, dt, wo, beta)

if(all == 1): #os comandos comentados são para exibição do gráfico usando a biblioteca pyplotly e outras funções mathplotlib

    # colors = ["blue", "green", "red", "orange"]
    # osc_cases = ["Sem Amortecimento", "Subamortecido", "Crítico", "Superamortecido"]

    # fig = go.Figure()
    #plt.style.use('dark_background')

######################################### Mude o valor de beta a partir daqui!. Atenção para as condições de cada caso

    print('Variável de controle na linha 7 ligada: \n' 
          'Executando todos os casos e plotando o gráfico de Posição x Tempo')

    beta = 0.0
    w = wo
    tempo, posicao, velocidade, ener_cin, ener_pot, energia = mhs_sa(xo, vo, m, k, npt, dt, wo, beta)

    plt.plot(tempo, posicao, label=rf'Sem Amortecimento, $\beta$ = {beta:.2f}', marker='o', color="blue")

    # fig.add_trace(go.Scatter(x=tempo, y=posicao,
    #                      mode='lines',
    #                      name=osc_cases[0],
    #                      line=dict(color=colors[0])))
    
    beta = wo - 0.9 * wo   
    w = math.sqrt(wo**2 - beta**2)
    tempo, posicao, velocidade, ener_cin, ener_pot, energia  = mhs_sub(xo, vo, m, k, npt, dt, w, beta)

    plt.plot(tempo, posicao, label=rf'Subamortecido, $\beta$ = {beta:.2f}', marker='^', color="green")

    # fig.add_trace(go.Scatter(x=tempo, y=posicao,
    #                      mode='lines',
    #                      name=osc_cases[1],
    #                      line=dict(color=colors[1])))

    beta = wo
    tempo, posicao, velocidade, ener_cin, ener_pot, energia = mhs_cri(xo, vo, m, k, npt, dt, wo, beta)

    plt.plot(tempo, posicao, label=rf'Crítico, $\beta = \omega$ = {beta:.2f}', marker='v', color="red")
   
    # fig.add_trace(go.Scatter(x=tempo, y=posicao,
    #                      mode='lines',
    #                      name=osc_cases[2],
    #                      line=dict(color=colors[2])))

    beta = wo + 0.9 * wo
    tempo, posicao, velocidade, ener_cin, ener_pot, energia = mhs_sup(xo, vo, m, k, npt, dt, wo, beta)
    
    # fig.add_trace(go.Scatter(x=tempo, y=posicao,
    #                      mode='lines',
    #                      name=osc_cases[3],
    #                      line=dict(color=colors[3])))
    # fig.update_layout(
    # title="Oscilações em Diferentes Casos",
    # xaxis_title="Tempo (s)",
    # yaxis_title="Amplitude (m)",
    # template="plotly_dark"
    # )

    # fig.show()

    plt.plot(tempo, posicao, label=rf'Superamortecido, $\beta$ = {beta:.2f}', marker='s', color="purple")

    plt.axhline(y=-0.000, color='black', linestyle='-')
    plt.xlabel('Tempo(s)', fontsize=14)
    plt.ylabel('Amplitude (m)', fontsize=14)
    plt.xlim(min(tempo), max(tempo))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.title('Movimento Harmônico Simples e Amortecido')
    plt.legend()
    plt.grid(False)
    plt.show()
# Salvando o gráfico
     #plt.savefig('grafico_mhs.png')

    # with open('saida_mhs.dat', 'w') as f:
    #     f.write("Tempo\tSem Amortecimento\tSubamortecido\tCrítico\tSuperamortecido\n")
    #     for i in range(npt):
    #         f.write(f"{tempo[i]:.4f}\t{posicao[i]:.4f}\t{posicao[i]:.4f}\t{posicao[i]:.4f}\t{posicao[i]:.4f}\n")

elif(all == 0):
    
    if (gr == 0):
                     
        plt.figure(figsize=(10, 6))
        plt.plot(tempo, posicao, label= 'Posição', marker='o', color="blue")
        plt.xlabel('Tempo(s)', fontsize=14)
        plt.ylabel('Amplitude (m)', labelpad=30.0, fontsize=14)
        plt.xlim(min(tempo)-0.01, max(tempo)+0.03)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.title(rf'(Posição x Tempo) para $\beta$ = {beta:.2f}', fontsize = 14)

        amp_max = max(posicao)
        amp_min = min(posicao)
        origem = 0
    
        plt.axhline(y=-0.000, color='black', linestyle='-')
        #plt.text(x=-1.5, y=-0.01, s=r"$O$", color='black', fontsize=18, va="center", ha="right")

        #plt.axhline(y=amp_max, color='red', label=f'Amplitude: {amp_max}')
        plt.text(x=-1.5, y=amp_max - 0.01, s=r"$+A$", color='black', fontsize=18, va="center", ha="right")

        #plt.axhline(y=amp_min, color='red', linestyle='--', label=f'Amplitude: {amp_min}')
        plt.text(x=-1.7, y=amp_min - 0.01, s=r"$-A$", color="black", fontsize=18, va="center", ha="right")

        plt.grid(False)
        plt.legend()
        plt.savefig('grafico_mhs_posicao.png')
        plt.show()        
    
    elif gr == 1:

        plt.figure(figsize=(10, 6))
        plt.plot(tempo, velocidade, label= 'Velocidade', marker='o', color="blue")
        plt.xlabel('Tempo(s)', fontsize=14)
        plt.ylabel('Velocidade (m/s)', labelpad=30.0, fontsize=14)
        plt.xlim(min(tempo)-0.01, max(tempo)+0.03)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.title(rf'(Velocidade x Tempo) para $\beta$ = {beta:.2f}', fontsize = 14)

        amp_max = max(velocidade)
        amp_min = min(velocidade)
        origem = 0
    
        plt.axhline(y=-0.000, color='black', linestyle='-')
        #plt.text(x=-1.5, y=-0.01, s=r"$O$", color='black', fontsize=18, va="center", ha="right")

        #plt.axhline(y=amp_max, color='red', label=f'Amplitude: {amp_max}')
        plt.text(x=-1.5, y=amp_max - 0.01, s=r"$+\omega A$", color='black', fontsize=18, va="center", ha="right")

        #plt.axhline(y=amp_min, color='red', linestyle='--', label=f'Amplitude: {amp_min}')
        plt.text(x=-1.7, y=amp_min - 0.01, s=r"$-\omega A$", color="black", fontsize=18, va="center", ha="right")

        plt.grid(False)
        plt.legend()
        plt.savefig('grafico_mhs_velocidade.png')
        plt.show()

    elif gr == 2:

        plt.figure(figsize=(10, 6))
        plt.plot(tempo, ener_cin, label= 'Energia Cinética', marker='o', color="green")
        plt.plot(tempo, ener_pot, label= 'Energia Potencial', marker='o', color="red")
        plt.xlabel('Tempo(s)', fontsize=14)
        plt.ylabel('Energia Total (J)', labelpad=30.0, fontsize=14)
        plt.xlim(min(tempo)-0.01, max(tempo)+0.03)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.title(rf'(Energia x Tempo) para $\beta$ = {beta:.2f}', fontsize = 14)

        amp_max = max(energia)
        amp_min = min(ener_cin)
        origem = 0
    
        plt.axhline(y=-0.000, color='black', linestyle='-')
        #plt.text(x=-1.5, y=-0.01, s=r"$O$", color='black', fontsize=18, va="center", ha="right")

        #plt.axhline(y=amp_max, color='red', label=f'Amplitude: {amp_max}')
        plt.text(x=-1.5, y=amp_max - 0.01, s=r"$Max$", color='black', fontsize=18, va="center", ha="right")

        #plt.axhline(y=amp_min, color='red', linestyle='--', label=f'Amplitude: {amp_min}')
        plt.text(x=-1.7, y=amp_min - 0.01, s=r"$Min$", color="black", fontsize=18, va="center", ha="right")

        plt.grid(False)
        plt.legend()
        plt.savefig('grafico_mhs_energia.png')
        plt.show()

if all not in (0,1):
    print('Erro: Verifique seus dados de entrada')
