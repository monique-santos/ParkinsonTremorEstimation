#include "sinal.h" 
#include <Arduino.h>
#include <Wire.h>
#include <Adafruit_VL53L0X.h>

// controle de impressão
bool imprimindo = false;
unsigned long inicioTeste = 0;
const unsigned long TEMPO_TESTE = 30000;

// sensor
Adafruit_VL53L0X lox = Adafruit_VL53L0X();

// ponte h
const int PIN_IN1 = 25;
const int PIN_IN2 = 14;
const int PIN_ENA = 33;

// pwm
const int frequenciaPWM = 500;
const int resolucaoPWM = 11;
int pwmReal = 0;

// i2c
const int PIN_SDA = 21;
const int PIN_SCL = 22;

// posição
float posicaoAtual = 0.0;
float posicaoInicialZero = 0.0;
float posicaoFiltrada = 0.0;

// filtro
const float ALPHA_FILTRO = 0.3;

// pid
float Kp = 3.0;
float Ki = 0.1;
float Kd = 0.0;

float erroAnterior = 0.0;
float integral = 0.0;

const float LIMITE_INTEGRAL = 100.0;

// pwm útil do motor
const int PWM_MINIMO = 1000;
const int PWM_MAX = 2047;

// tolerância
const float ZONA_MORTA = 1.0;

// tempo do loop
unsigned long ultimoTempoLoop = 0;
const unsigned long INTERVALO_SISTEMA = 20;

// referência
int indiceReferencia = 0;
float referenciaAtual = 0.0;
int tam_referencia = TAMANHO_SINAL;

void setup() {

    Serial.begin(115200);

    Wire.begin(PIN_SDA, PIN_SCL);

    configurarPinos();
    configurarPWM_ESP32();

    inicializarLaser();

    calibrarPosicaoInicial();

    ultimoTempoLoop = millis();

    imprimindo = true;
    inicioTeste = millis();

    Serial.println("Tempo_ms,Posicao_mm,Controle,PWM_Real,Referencia_mm");
}

void loop() {

    unsigned long tempoAtual = millis();

    // loop fixo em 20 ms
    if (tempoAtual - ultimoTempoLoop >= INTERVALO_SISTEMA) {

        ultimoTempoLoop = tempoAtual;

        atualizarReferencia();

        processarPID(tempoAtual);
    }
}

void configurarPinos() {

    pinMode(PIN_IN1, OUTPUT);
    pinMode(PIN_IN2, OUTPUT);

    digitalWrite(PIN_IN1, LOW);
    digitalWrite(PIN_IN2, LOW);
}

void configurarPWM_ESP32() {

    ledcAttach(PIN_ENA, frequenciaPWM, resolucaoPWM);

    ledcWrite(PIN_ENA, 0);
}

void inicializarLaser() {

    delay(100);

    if (!lox.begin()) {

        while (1);
    }

    // leitura em 20 ms
    lox.setMeasurementTimingBudgetMicroSeconds(20000);
}

void calibrarPosicaoInicial() {

    pararMotor();

    float soma = 0;
    int amostras = 0;

    unsigned long start = millis();

    // média inicial
    while (millis() - start < 5000) {

        VL53L0X_RangingMeasurementData_t measure;

        lox.rangingTest(&measure, false);

        if (measure.RangeStatus != 4) {

            soma += measure.RangeMilliMeter;

            amostras++;
        }

        delay(20);
    }

    posicaoInicialZero =
        (amostras > 0) ?
        (soma / amostras) :
        0;

    posicaoFiltrada = 0;
}

void atualizarReferencia() {

    if (tam_referencia > 0) {

        referenciaAtual =
            sinal_referencia[indiceReferencia] * 1.0;

        indiceReferencia =
            (indiceReferencia + 1) % tam_referencia;
    }
}

void processarPID(unsigned long tempoAtual) {

    float dt = 0.02;

    VL53L0X_RangingMeasurementData_t measure;

    lox.rangingTest(&measure, false);

    if (measure.RangeStatus != 4) {

        // posição relativa ao zero
        float leituraBruta =
            measure.RangeMilliMeter - posicaoInicialZero;

        // filtro
        posicaoFiltrada =
            (ALPHA_FILTRO * leituraBruta) +
            ((1.0 - ALPHA_FILTRO) * posicaoFiltrada);

        posicaoAtual = posicaoFiltrada;

        // erro
        float erro =
            referenciaAtual - posicaoAtual;

        // pid
        float controle =
            calcularPID(erro, dt);

        // motor
        acionarMotor(controle, erro);

        // serial
        if (imprimindo) {

            unsigned long tempoRelativo =
                millis() - inicioTeste;

            Serial.print(tempoRelativo);
            Serial.print(",");

            Serial.print(posicaoAtual);
            Serial.print(",");

            Serial.print(controle);
            Serial.print(",");

            Serial.print(pwmReal);
            Serial.print(",");

            Serial.println(referenciaAtual);

            if (tempoRelativo >= TEMPO_TESTE) {

                imprimindo = false;

                pararMotor();
            }
        }

    } else {

        pararMotor();
    }
}

float calcularPID(float erro, float dt) {

    // proporcional
    float P = Kp * erro;

    // integral
    integral += erro * dt;

    integral =
        constrain(
            integral,
            -LIMITE_INTEGRAL,
            LIMITE_INTEGRAL
        );

    float I = Ki * integral;

    // derivativo
    float D =
        (erro - erroAnterior) / dt;

    erroAnterior = erro;

    return P + I + D;
}

void acionarMotor(float sinalControle, float erro) {

    // zona morta
    if (abs(erro) < ZONA_MORTA) {

        pararMotor();

        return;
    }

    // módulo do controle
    float absControle =
        abs(sinalControle);

    absControle =
        constrain(
            absControle,
            0.0,
            100.0
        );

    // conversão pra pwm
    int pwmValue =
        PWM_MINIMO +
        (int)(
            (absControle / 100.0) *
            (PWM_MAX - PWM_MINIMO)
        );

    pwmValue =
        constrain(
            pwmValue,
            PWM_MINIMO,
            PWM_MAX
        );

    pwmReal = pwmValue;

    // direção
    if (sinalControle > 0) {

        moverFrente(pwmValue);

    } else {

        moverTras(pwmValue);
    }
}

void moverFrente(int velocidade) {

    digitalWrite(PIN_IN1, HIGH);
    digitalWrite(PIN_IN2, LOW);

    ledcWrite(PIN_ENA, velocidade);
}

void moverTras(int velocidade) {

    digitalWrite(PIN_IN1, LOW);
    digitalWrite(PIN_IN2, HIGH);

    ledcWrite(PIN_ENA, velocidade);
}

void pararMotor() {

    digitalWrite(PIN_IN1, LOW);
    digitalWrite(PIN_IN2, LOW);

    ledcWrite(PIN_ENA, 0);

    pwmReal = 0;
}