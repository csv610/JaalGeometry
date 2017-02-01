/****************************************************************************
** Meta object code from reading C++ file 'MeshWavefrontsDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.5.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "MeshWavefrontsDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MeshWavefrontsDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.5.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JMeshWavefrontsDialog_t {
    QByteArrayData data[9];
    char stringdata0[96];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JMeshWavefrontsDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JMeshWavefrontsDialog_t qt_meta_stringdata_JMeshWavefrontsDialog = {
    {
        QT_MOC_LITERAL(0, 0, 21), // "JMeshWavefrontsDialog"
        QT_MOC_LITERAL(1, 22, 8), // "nextWave"
        QT_MOC_LITERAL(2, 31, 0), // ""
        QT_MOC_LITERAL(3, 32, 13), // "waveAnimation"
        QT_MOC_LITERAL(4, 46, 11), // "closeDialog"
        QT_MOC_LITERAL(5, 58, 10), // "setNewWave"
        QT_MOC_LITERAL(6, 69, 13), // "keyPressEvent"
        QT_MOC_LITERAL(7, 83, 10), // "QKeyEvent*"
        QT_MOC_LITERAL(8, 94, 1) // "e"

    },
    "JMeshWavefrontsDialog\0nextWave\0\0"
    "waveAnimation\0closeDialog\0setNewWave\0"
    "keyPressEvent\0QKeyEvent*\0e"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JMeshWavefrontsDialog[] = {

// content:
    7,       // revision
    0,       // classname
    0,    0, // classinfo
    5,   14, // methods
    0,    0, // properties
    0,    0, // enums/sets
    0,    0, // constructors
    0,       // flags
    0,       // signalCount

// slots: name, argc, parameters, tag, flags
    1,    0,   39,    2, 0x08 /* Private */,
    3,    0,   40,    2, 0x08 /* Private */,
    4,    0,   41,    2, 0x08 /* Private */,
    5,    0,   42,    2, 0x08 /* Private */,
    6,    1,   43,    2, 0x08 /* Private */,

// slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 7,    8,

    0        // eod
};

void JMeshWavefrontsDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JMeshWavefrontsDialog *_t = static_cast<JMeshWavefrontsDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0:
            _t->nextWave();
            break;
        case 1:
            _t->waveAnimation();
            break;
        case 2:
            _t->closeDialog();
            break;
        case 3:
            _t->setNewWave();
            break;
        case 4:
            _t->keyPressEvent((*reinterpret_cast< QKeyEvent*(*)>(_a[1])));
            break;
        default:
            ;
        }
    }
}

const QMetaObject JMeshWavefrontsDialog::staticMetaObject = {
    {   &QDialog::staticMetaObject, qt_meta_stringdata_JMeshWavefrontsDialog.data,
        qt_meta_data_JMeshWavefrontsDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR
    }
};


const QMetaObject *JMeshWavefrontsDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JMeshWavefrontsDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JMeshWavefrontsDialog.stringdata0))
        return static_cast<void*>(const_cast< JMeshWavefrontsDialog*>(this));
    if (!strcmp(_clname, "Ui::MeshWavefrontsDialog"))
        return static_cast< Ui::MeshWavefrontsDialog*>(const_cast< JMeshWavefrontsDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JMeshWavefrontsDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 5)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 5;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 5)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 5;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
