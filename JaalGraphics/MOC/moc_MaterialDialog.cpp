/****************************************************************************
** Meta object code from reading C++ file 'MaterialDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.5.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "MaterialDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MaterialDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.5.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JMaterialDialog_t {
    QByteArrayData data[7];
    char stringdata0[78];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JMaterialDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JMaterialDialog_t qt_meta_stringdata_JMaterialDialog = {
    {
        QT_MOC_LITERAL(0, 0, 15), // "JMaterialDialog"
        QT_MOC_LITERAL(1, 16, 10), // "setAmbient"
        QT_MOC_LITERAL(2, 27, 0), // ""
        QT_MOC_LITERAL(3, 28, 10), // "setDiffuse"
        QT_MOC_LITERAL(4, 39, 11), // "setEmission"
        QT_MOC_LITERAL(5, 51, 11), // "setSpecular"
        QT_MOC_LITERAL(6, 63, 14) // "setStdMaterial"

    },
    "JMaterialDialog\0setAmbient\0\0setDiffuse\0"
    "setEmission\0setSpecular\0setStdMaterial"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JMaterialDialog[] = {

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
    6,    0,   43,    2, 0x08 /* Private */,

// slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

    0        // eod
};

void JMaterialDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JMaterialDialog *_t = static_cast<JMaterialDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0:
            _t->setAmbient();
            break;
        case 1:
            _t->setDiffuse();
            break;
        case 2:
            _t->setEmission();
            break;
        case 3:
            _t->setSpecular();
            break;
        case 4:
            _t->setStdMaterial();
            break;
        default:
            ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject JMaterialDialog::staticMetaObject = {
    {   &QDialog::staticMetaObject, qt_meta_stringdata_JMaterialDialog.data,
        qt_meta_data_JMaterialDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR
    }
};


const QMetaObject *JMaterialDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JMaterialDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JMaterialDialog.stringdata0))
        return static_cast<void*>(const_cast< JMaterialDialog*>(this));
    if (!strcmp(_clname, "Ui::MaterialDialog"))
        return static_cast< Ui::MaterialDialog*>(const_cast< JMaterialDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JMaterialDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
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
